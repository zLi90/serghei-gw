/* -*- mode: c++; c-default-style: "linux" -*- */

/*******************************************************************************
 * @file    MeteoInit.h
 * @brief   Meteorological data reader and reference evapotranspiration calculator
 *          for the WOFOST crop growth model coupled with SERGHEI.
 *
 * @details This module performs two primary functions:
 *
 *   1. **Data ingestion**: Reads a daily meteorological forcing file that
 *      contains columns for time, year, day-of-year (DOY), solar irradiation,
 *      minimum temperature, maximum temperature, vapour pressure, wind speed,
 *      and rainfall. It also extracts site metadata (latitude, longitude,
 *      elevation, and optional Angstrom coefficients) from header comments.
 *
 *   2. **Reference ET computation**: For each daily record, computes three
 *      reference evapotranspiration rates using established methods:
 *      - E0  : Open water evaporation (Penman 1948)
 *      - ES0 : Bare soil evaporation (Penman 1948)
 *      - ET0 : Reference crop evapotranspiration (Penman-Monteith / FAO-56 style)
 *      All three are stored in units of [cm/day].
 *
 *          The meteorological file format is expected to be:
 *            - Header lines starting with "//" that may contain site info:
 *                // LAT 32.79
 *                // LON 119.45
 *                // ELEV 7
 *                // ANGSTROMA 0.18
 *                // ANGSTROMB 0.50
 *            - Data lines with 9 or more columns:
 *                Time  Year  DOY  Irrad  Tmin  Tmax  Vap  Wind  Rain
 *                (1.0  2000  1    5.6e6  -2.0  5.0   6.5  3.2   0.0)
 *
 *          Column units:
 *            - Irrad: Global daily solar radiation [J/m^2/d]
 *            - Tmin, Tmax: Air temperature [degrees Celsius]
 *            - Vap: Vapour pressure [hPa]
 *            - Wind: Wind speed at 2 m height [m/s]
 *            - Rain: Daily precipitation [mm/day or equivalent]
 *
 * @see     MeteoState.h  for the state container that stores all meteorological arrays
 ******************************************************************************/

#ifndef _METEO_INIT_H_
#define _METEO_INIT_H_

#include "../define.h"
#include "../Parallel.h"
#include "MeteoState.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <cmath>

/**
 * @class MeteoInit
 * @brief Reads meteorological forcing data and computes reference evapotranspiration.
 *
 * This class provides functionality to:
 *   - Parse a text-based meteorological file with header metadata
 *   - Extract site location parameters (latitude, longitude, elevation)
 *   - Optionally parse Angstrom A and B coefficients for radiation conversion
 *   - Compute daily E0, ES0, and ET0 using the Penman (1948) and
 *     Penman-Monteith methods
 *   - Store all data into Kokkos device views via MeteoState
 */
class MeteoInit
{

public:
    /* ======================================================================= *
     * Primary Data Reader
     * ======================================================================= */

    /**
     * @brief Read meteorological data from file and populate MeteoState.
     *
     * This method opens the specified file, parses location metadata from
     * comment headers, reads all daily data records, computes reference
     * evapotranspiration for each day, and transfers everything to Kokkos
     * device memory.
     *
     * @param[in]  fNameIn  Path to the meteorological input file
     * @param[out] ms       MeteoState object to populate with data and ET values
     * @param[in]  par      Parallel context (used for master-process-only I/O)
     * @return 1 on success, 0 on failure
     */
    int readMeteoFile(std::string fNameIn, MeteoState &ms, Parallel &par)
    {
        std::ifstream fInStream(fNameIn);

        if (!fInStream.is_open())
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Unable to open meteo file: " << fNameIn << "\n";
            }
            return 0;
        }

        std::string line;

        /* Temporary storage on Host using std::vector.
         * Each inner vector holds one row of numeric data columns. */
        std::vector<std::vector<real>> raw_data;

        /* ----------------------------------------------------------------- *
         * Default location parameters.
         * The sentinel value -999.0 indicates "not yet parsed".
         * These will be overwritten if the file contains LAT/LON/ELEV
         * header comments; otherwise an error is raised below.
         * ----------------------------------------------------------------- */
        ms.lat = -999.0;    /* Latitude [decimal degrees], range: -90 to +90 */
        ms.lon = -999.0;    /* Longitude [decimal degrees], range: -180 to +180 */
        ms.elev = -999.0;   /* Elevation above sea level [m], range: typically -100 to +5000 */

        /* ----------------------------------------------------------------- *
         * Default Angstrom coefficients for estimating sunshine duration
         * from measured global radiation.
         *   angstA: Angstrom A coefficient [dimensionless], typically 0.10 - 0.30
         *   angstB: Angstrom B coefficient [dimensionless], typically 0.40 - 0.70
         * These defaults (0.18, 0.50) correspond to clear-sky temperate conditions.
         * ----------------------------------------------------------------- */
        ms.angstA = 0.18;
        ms.angstB = 0.50;

        /* Store all file lines for two-pass parsing (location info + data). */
        std::vector<std::string> file_lines;

        while (std::getline(fInStream, line))
        {
            file_lines.push_back(line);

            /* Trim leading whitespace to detect line type. */
            size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos)
                continue; /* Skip blank lines */

            /* Attempt to parse location metadata from this line.
             * Works on both comment lines (// LAT ...) and plain lines. */
            parseLocationInfo(line, ms);

            /* Skip full-line comments (starting with //) for data parsing. */
            if (line[first] == '/' && line.size() > first + 1 && line[first + 1] == '/')
            {
                continue;
            }

            /* Data lines: must start with a digit or a minus sign. */
            if (isdigit(line[first]) || line[first] == '-')
            {
                std::stringstream ss(line);
                std::vector<real> row_vals;
                real val;
                while (ss >> val)
                {
                    row_vals.push_back(val);
                }

                /* Expected columns: Time, Year, DOY, Irrad, Tmin, Tmax, Vap, Wind, Rain
                 * Total 9 columns minimum required.                  0     1     2    3      4     5     6    7     8   */
                if (row_vals.size() >= 9)
                {
                    raw_data.push_back(row_vals);
                }
            }
        }
        fInStream.close();

        /* =================================================================== *
         * Validation: Verify that location parameters were successfully parsed
         * =================================================================== */
        bool location_error = false;
        std::string error_msg = "";

        if (ms.lat == -999.0)
        {
            location_error = true;
            error_msg += "Latitude (lat) was not correctly read from meteo file or is still -999.\n";
        }

        if (ms.lon == -999.0)
        {
            location_error = true;
            error_msg += "Longitude (lon) was not correctly read from meteo file or is still -999.\n";
        }

        if (ms.elev == -999.0)
        {
            location_error = true;
            error_msg += "Elevation (elev) was not correctly read from meteo file or is still -999.\n";
        }

        if (location_error)
        {
            if (par.masterproc)
            {
                std::cerr << RERROR "Failed to read location parameters from meteo file: " << fNameIn << "\n";
                std::cerr << RERROR << error_msg;
                std::cerr << RERROR "Please ensure the meteo file contains valid LAT, LON, and ELEV values (e.g., '// LAT 32.79').\n";
            }
            return 0;
        }

        if (par.masterproc)
        {
            std::cout << GOK "Location parameters verified - Lat=" << ms.lat
                      << ", Lon=" << ms.lon << ", Elev=" << ms.elev << std::endl;
        }

        /* =================================================================== *
         * Data validation: Ensure at least one valid record was found
         * =================================================================== */
        size_t n = raw_data.size();
        if (n == 0)
        {
            if (par.masterproc)
                std::cerr << RERROR "No valid meteo records found.\n";
            return 0;
        }

        /* =================================================================== *
         * Extract simulation start date from the first data record
         *
         * Column 1 = Year (e.g., 2000)
         * Column 2 = Day of Year / DOY (1-366)
         * =================================================================== */
        if (raw_data[0].size() >= 3)
        {
            ms.start_year = (int)raw_data[0][1];
            ms.start_doy = (int)raw_data[0][2];

            if (par.masterproc)
            {
                std::cout << GOK "Meteo Start Date Detected: Year " << ms.start_year
                          << ", DOY " << ms.start_doy << std::endl;
            }
        }
        else
        {
            /* Fallback to a safe default if columns are missing. */
            ms.start_year = 2000;
            ms.start_doy = 1;
        }

        /* =================================================================== *
         * Allocate Kokkos Views on Device and populate with parsed data
         * =================================================================== */

        /* Step 1: Allocate device-side Kokkos Views of size n. */
        ms.allocate(n);

        /* Step 2: Create host-side mirror views for CPU-accessible writes. */
        auto h_time = Kokkos::create_mirror_view(ms.time);
        auto h_irrad = Kokkos::create_mirror_view(ms.irrad);
        auto h_tmin = Kokkos::create_mirror_view(ms.tmin);
        auto h_tmax = Kokkos::create_mirror_view(ms.tmax);
        auto h_vap = Kokkos::create_mirror_view(ms.vap);
        auto h_wind = Kokkos::create_mirror_view(ms.wind);
        auto h_rain = Kokkos::create_mirror_view(ms.rain);

        /* Mirror views for derived ET variables (initialized to 0 by Kokkos). */
        auto h_e0 = Kokkos::create_mirror_view(ms.e0);
        auto h_es0 = Kokkos::create_mirror_view(ms.es0);
        auto h_et0 = Kokkos::create_mirror_view(ms.et0);

        /* Step 3: Fill host mirrors with parsed data and compute reference ET. */
        for (size_t i = 0; i < n; ++i)
        {
            /* ----------------------------------------------------------------- *
             * Column mapping from the raw data file:
             *   Index 0: Time          (sequential day counter or simulation time)
             *   Index 1: Year          (calendar year, e.g., 2000)
             *   Index 2: DOY           (day of year, 1-366)
             *   Index 3: Irrad         (global solar radiation [J/m^2/d])
             *   Index 4: Tmin          (daily minimum temperature [deg C])
             *   Index 5: Tmax          (daily maximum temperature [deg C])
             *   Index 6: Vap           (vapour pressure [hPa])
             *   Index 7: Wind          (wind speed at 2 m height [m/s])
             *   Index 8: Rain          (daily precipitation [mm/day])
             * ----------------------------------------------------------------- */
            h_time(i) = raw_data[i][0];
            h_irrad(i) = raw_data[i][3];
            h_tmin(i) = raw_data[i][4];
            h_tmax(i) = raw_data[i][5];
            h_vap(i) = raw_data[i][6];
            h_wind(i) = raw_data[i][7];
            h_rain(i) = raw_data[i][8];

            /* ----------------------------------------------------------------- *
             * Prepare inputs for the ET calculation.
             *
             * Irradiation:  [J/m^2/d] - total daily shortwave radiation
             *   Typical range: 0 (polar night) to 3.5e7 J/m^2/d (tropical summer)
             *
             * Vapour pressure: [hPa] - atmospheric water vapour pressure
             *   Typical range: 1 hPa (cold/dry) to 40 hPa (hot/humid)
             *
             * Wind speed: [m/s] - measured at 2 m above ground
             *   Typical range: 0.5 to 10 m/s
             * ----------------------------------------------------------------- */
            real TMIN = raw_data[i][4];
            real TMAX = raw_data[i][5];
            real IRRAD = raw_data[i][3]; /* [J/m^2/d] */
            real VAP = raw_data[i][6];   /* [hPa] */
            real WIND = raw_data[i][7];  /* [m/s] */

            /* ----------------------------------------------------------------- *
             * Convert Year and DOY to a std::tm struct for astronomical
             * calculations (solar declination, day length, etc.).
             *
             * The approach constructs Jan 1 of the current year, then
             * advances by (DOY - 1) days using mktime normalization.
             * ----------------------------------------------------------------- */
            int current_year = (int)raw_data[i][1];
            int current_doy = (int)raw_data[i][2];

            std::tm date_tm = {};
            date_tm.tm_year = current_year - 1900;  /* tm_year is years since 1900 */
            date_tm.tm_mday = 1;                     /* Start from January 1 */
            date_tm.tm_mon = 0;                      /* January (0-based) */
            std::mktime(&date_tm);                   /* Normalize the tm struct */

            /* Advance by (DOY - 1) days since tm_yday is 0-based.
             * mktime handles month/day overflow automatically. */
            date_tm.tm_mday += (current_doy - 1);
            std::mktime(&date_tm);                   /* Re-normalize */

            /* ----------------------------------------------------------------- *
             * Compute reference evapotranspiration using combined Penman (1948)
             * for E0/ES0 and Penman-Monteith for ET0.
             *
             * Output units from the calculation: [mm/day]
             * Stored units after conversion:     [cm/day]
             * Conversion factor: 0.1 (1 mm = 0.1 cm)
             *
             * Typical ranges for reference ET0 (temperate climate):
             *   Winter: 0.05 - 0.20 cm/day  (0.5 - 2.0 mm/day)
             *   Summer: 0.30 - 0.70 cm/day  (3.0 - 7.0 mm/day)
             * ----------------------------------------------------------------- */
            real e0, es0, et0;
            calculateET(date_tm, ms.lat, ms.elev, TMIN, TMAX, IRRAD, VAP, WIND,
                        ms.angstA, ms.angstB, e0, es0, et0);

            /* Convert from mm/day to cm/day and store directly on device views. */
            ms.e0(i) = e0 * 0.1;    /* Open water evaporation [cm/day] */
            ms.es0(i) = es0 * 0.1;  /* Bare soil evaporation [cm/day] */
            ms.et0(i) = et0 * 0.1;  /* Reference crop evapotranspiration [cm/day] */

#if DEBUG_CROP_GROWTH_MODEL
            printf("MeteoInit: i=%zu, DOY=%d, E0=%.6f cm/day, ES0=%.6f cm/day, ET0=%.6f cm/day\n",
                   i, date_tm.tm_mday, h_e0(i), h_es0(i), ms.et0(i));
#endif
        }

        /* Step 4: Deep copy host mirror data to Kokkos device views. */
        Kokkos::deep_copy(ms.time, h_time);
        Kokkos::deep_copy(ms.irrad, h_irrad);
        Kokkos::deep_copy(ms.tmin, h_tmin);
        Kokkos::deep_copy(ms.tmax, h_tmax);
        Kokkos::deep_copy(ms.vap, h_vap);
        Kokkos::deep_copy(ms.wind, h_wind);
        Kokkos::deep_copy(ms.rain, h_rain);

#if DEBUG_CROP_GROWTH_MODEL
        /* Debug: Copy irradiation back to host and print all values. */
        {
            auto h_irrad_out = Kokkos::create_mirror_view(ms.irrad);
            Kokkos::deep_copy(h_irrad_out, ms.irrad);
            for (size_t i = 0; i < n; ++i)
            {
                printf("ms.irrad[%zu] = %.6f cm/day\n", i, static_cast<double>(h_irrad_out(i)));
            }
        }
#endif

        if (par.masterproc)
        {
            std::cout << GOK "Crop Meteo data read from " << fNameIn << "\n";
            std::cout << GOK "Crop Meteo data read: " << n << " records." << std::endl;
        }

        return 1;
    }

private:
    /* ======================================================================= *
     * Location Metadata Parser
     * ======================================================================= */

    /**
     * @brief Parse site location and radiation parameters from comment lines.
     *
     * Scans a single text line for keywords LAT, LON, ELEV, ANGSTROMA, ANGSTROMB
     * (case-insensitive) and extracts the associated numeric value. Handles
     * optional colon separators (e.g., "LAT: 32.79" or "LAT 32.79").
     *
     * @param[in]     line  Input text line to parse
     * @param[in,out] ms    MeteoState to store parsed location parameters
     */
    void parseLocationInfo(const std::string &line, MeteoState &ms)
    {
        std::string trimmed = line;

        /* Remove leading "//" comment marker if present. */
        size_t pos = trimmed.find("//");
        if (pos != std::string::npos)
        {
            trimmed = trimmed.substr(pos + 2);
        }

        /* Trim leading and trailing whitespace. */
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

        /* Convert to uppercase for case-insensitive keyword matching. */
        std::string upperLine = trimmed;
        std::transform(upperLine.begin(), upperLine.end(), upperLine.begin(), ::toupper);

        /* --- Parse latitude [decimal degrees] --- */
        if (upperLine.find("LAT") != std::string::npos)
        {
            size_t lat_pos = upperLine.find("LAT");
            std::string substr = trimmed.substr(lat_pos + 3);
            /* Skip optional colon and whitespace after keyword. */
            size_t first = substr.find_first_not_of(" \t:");
            if (first != std::string::npos)
            {
                substr = substr.substr(first);
                std::stringstream ss(substr);
                double val;
                if (ss >> val)
                {
                    ms.lat = val;
                }
            }
        }

        /* --- Parse longitude [decimal degrees] --- */
        if (upperLine.find("LON") != std::string::npos)
        {
            size_t lon_pos = upperLine.find("LON");
            std::string substr = trimmed.substr(lon_pos + 3);
            size_t first = substr.find_first_not_of(" \t:");
            if (first != std::string::npos)
            {
                substr = substr.substr(first);
                std::stringstream ss(substr);
                double val;
                if (ss >> val)
                {
                    ms.lon = val;
                }
            }
        }

        /* --- Parse elevation [m above sea level] --- */
        if (upperLine.find("ELEV") != std::string::npos)
        {
            size_t elev_pos = upperLine.find("ELEV");
            std::string substr = trimmed.substr(elev_pos + 4);
            size_t first = substr.find_first_not_of(" \t:");
            if (first != std::string::npos)
            {
                substr = substr.substr(first);
                std::stringstream ss(substr);
                double val;
                if (ss >> val)
                {
                    ms.elev = val;
                }
            }
        }

        /* --- Parse Angstrom A coefficient [dimensionless, typical: 0.10 - 0.30] --- */
        if (upperLine.find("ANGSTROMA") != std::string::npos)
        {
            size_t anga_pos = upperLine.find("ANGSTROMA");
            std::string substr = trimmed.substr(anga_pos + 9);
            size_t first = substr.find_first_not_of(" \t:");
            if (first != std::string::npos)
            {
                substr = substr.substr(first);
                std::stringstream ss(substr);
                double val;
                if (ss >> val)
                {
                    ms.angstA = val;
                }
            }
        }

        /* --- Parse Angstrom B coefficient [dimensionless, typical: 0.40 - 0.70] --- */
        if (upperLine.find("ANGSTROMB") != std::string::npos)
        {
            size_t angb_pos = upperLine.find("ANGSTROMB");
            std::string substr = trimmed.substr(angb_pos + 9);
            size_t first = substr.find_first_not_of(" \t:");
            if (first != std::string::npos)
            {
                substr = substr.substr(first);
                std::stringstream ss(substr);
                double val;
                if (ss >> val)
                {
                    ms.angstB = val;
                }
            }
        }
        printf("MeteoInit: Parsed Location Info - Lat=%.2f, Lon=%.2f, Elev=%.1f, AngstromA=%.3f, AngstromB=%.3f\n",
               ms.lat, ms.lon, ms.elev, ms.angstA, ms.angstB);
    }

    /* ======================================================================= *
     * Unit Conversion and Thermodynamic Helper Functions
     * ======================================================================= */

    /**
     * @brief Convert pressure from hectopascals to kilopascals.
     * @param hPa  Pressure in hectopascals [hPa]
     * @return Pressure in kilopascals [kPa]
     */
    real hPa2kPa(real hPa) const
    {
        return hPa / 10.0;
    }

    /**
     * @brief Convert temperature from degrees Celsius to Kelvin.
     * @param celsius  Temperature in degrees Celsius
     * @return Temperature in Kelvin [K]
     */
    real Celsius2Kelvin(real celsius) const
    {
        return celsius + 273.16;
    }

    /**
     * @brief Compute saturation vapour pressure using the Tetens formula.
     *
     * Uses the empirical relation:
     *   e_s(T) = 0.6108 * exp(17.27 * T / (237.3 + T))
     *
     * @param temp  Air temperature [degrees Celsius]
     * @return Saturated vapour pressure [kPa]
     *
     * Typical range: 0.61 kPa at 0 deg C to 7.38 kPa at 40 deg C
     */
    real SatVapourPressure(real temp) const
    {
        return 0.6108 * exp((17.27 * temp) / (237.3 + temp));
    }

    /**
     * @brief Clamp a value to the range [vmin, vmax].
     * @param vmin  Minimum allowed value
     * @param vmax  Maximum allowed value
     * @param v     Input value
     * @return Clamped value
     */
    real limit(real vmin, real vmax, real v) const
    {
        if (v < vmin)
            return vmin;
        if (v > vmax)
            return vmax;
        return v;
    }

    /* ======================================================================= *
     * Astronomical Calculations (Solar Geometry)
     * ======================================================================= */

    /**
     * @brief Result structure for astronomical calculations.
     *
     * Contains solar geometry parameters computed for a given day and latitude.
     */
    struct AstroResults
    {
        real DAYL;    /* Day length [hours], range: 0 - 24 */
        real DAYLP;   /* Day length with civil twilight correction (base angle = -4 deg) [hours] */
        real SINLD;   /* sin(latitude) * sin(declination) [dimensionless] */
        real COSLD;   /* cos(latitude) * cos(declination) [dimensionless] */
        real DIFPP;   /* Diffuse radiation contribution parameter [J/m^2/s or similar] */
        real ATMTR;   /* Atmospheric transmission coefficient [dimensionless], range: 0 - 1+ */
        real DSINBE;  /* Integrated sin(β) over the day, used for radiation integration [s] */
        real ANGOT;   /* Angot extraterrestrial radiation (top of atmosphere) [J/m^2/d] */
    };

    /**
     * @brief Compute solar geometry parameters for a given day and latitude.
     *
     * Calculates day length, solar declination, extraterrestrial radiation,
     * atmospheric transmission, and diffuse radiation fraction using standard
     * astronomical equations.
     *
     * @param day_of_year  Day of year (1-366)
     * @param latitude     Site latitude [decimal degrees, positive = North]
     * @param radiation    Measured global solar radiation at surface [J/m^2/d]
     * @return AstroResults structure with all computed parameters
     */
    AstroResults astro(int day_of_year, real latitude, real radiation) const
    {
        AstroResults result = {};

        /* Reject invalid latitude values. */
        if (abs(latitude) > 90.0)
        {
            return result;
        }

        real RAD = M_PI / 180.0;  /* Degrees to radians conversion factor */
        real LAT = latitude;

        /* --- Solar declination ---
         * The declination angle varies sinusoidally over the year.
         * At DOY=172 (summer solstice in NH), declination ~ +23.45 deg.
         * At DOY=355 (winter solstice in NH), declination ~ -23.45 deg. */
        real DEC = -asin(sin(23.45 * RAD) * cos(2.0 * M_PI * (day_of_year + 10.0) / 365.0));

        /* --- Solar constant corrected for Earth-Sun distance ---
         * Base solar constant: 1370 W/m^2 (the Scafetta & Willson 2019 value is ~1361, but
         * 1370 is the WOFOST convention). Eccentricity correction applied. */
        real SC = 1370.0 * (1.0 + 0.033 * cos(2.0 * M_PI * day_of_year / 365.0));

        /* --- Intermediate variables for sunrise/sunset calculations --- */
        real SINLD = sin(RAD * LAT) * sin(DEC);
        real COSLD = cos(RAD * LAT) * cos(DEC);
        real AOB = SINLD / COSLD;

        /* --- Day length and integrated solar angles ---
         * DAYL: hours between sunrise and sunset (sun center at horizon, 0 deg).
         * DSINB: integrated sin(solar elevation) over the day [s].
         * DSINBE: same but weighted for diffuse radiation calculations [s]. */
        real DAYL = 0.0;
        real DSINB = 0.0;
        real DSINBE = 0.0;

        if (abs(AOB) <= 1.0)
        {
            /* Normal case: sunrise and sunset occur within the day. */
            DAYL = 12.0 * (1.0 + 2.0 * asin(AOB) / M_PI);
            DSINB = 3600.0 * (DAYL * SINLD + 24.0 * COSLD * sqrt(1.0 - AOB * AOB) / M_PI);
            DSINBE = 3600.0 * (DAYL * (SINLD + 0.4 * (SINLD * SINLD + COSLD * COSLD * 0.5)) +
                               12.0 * COSLD * (2.0 + 3.0 * 0.4 * SINLD) * sqrt(1.0 - AOB * AOB) / M_PI);
        }
        else
        {
            /* Polar cases: midnight sun or polar night. */
            if (AOB > 1.0)
                DAYL = 24.0;   /* Sun never sets */
            if (AOB < -1.0)
                DAYL = 0.0;    /* Sun never rises */
            DSINB = 3600.0 * (DAYL * SINLD);
            DSINBE = 3600.0 * (DAYL * (SINLD + 0.4 * (SINLD * SINLD + COSLD * COSLD * 0.5)));
        }

        /* --- Day length with civil twilight (sun at -4 deg below horizon) ---
         * This is used to extend the effective day length for certain
         * photoperiod-sensitive crop calculations. */
        real ANGLE = -4.0;
        real AOB_CORR = (-sin(ANGLE * RAD) + SINLD) / COSLD;
        real DAYLP = 0.0;

        if (abs(AOB_CORR) <= 1.0)
        {
            DAYLP = 12.0 * (1.0 + 2.0 * asin(AOB_CORR) / M_PI);
        }
        else if (AOB_CORR > 1.0)
        {
            DAYLP = 24.0;
        }
        else if (AOB_CORR < -1.0)
        {
            DAYLP = 0.0;
        }

        /* --- Extraterrestrial (Angot) radiation and atmospheric transmission ---
         * ANGOT: Total daily extraterrestrial radiation at the top of atmosphere [J/m^2/d]
         * ATMTR: Ratio of measured to extraterrestrial radiation [dimensionless]
         *   Typical clear sky: 0.7 - 0.8
         *   Overcast:          0.2 - 0.4 */
        real ANGOT = SC * DSINB;
        real ATMTR = 0.0;
        if (DAYL > 0.0)
        {
            ATMTR = radiation / ANGOT;
        }

        /* --- Estimate fraction of diffuse irradiation ---
         * Based on the atmospheric transmission ratio (ATMTR):
         *   High ATMTR (>0.75): Clear sky, low diffuse fraction (0.23)
         *   Moderate (0.35-0.75): Partial cloud, increasing diffuse
         *   Low ATMTR (0.07-0.35): Overcast, mostly diffuse
         *   Very low (<0.07): Nearly all diffuse (1.0) */
        real FRDIF = 0.0;
        if (ATMTR > 0.75)
        {
            FRDIF = 0.23;
        }
        else if (ATMTR <= 0.75 && ATMTR > 0.35)
        {
            FRDIF = 1.33 - 1.46 * ATMTR;
        }
        else if (ATMTR <= 0.35 && ATMTR > 0.07)
        {
            FRDIF = 1.0 - 2.3 * (ATMTR - 0.07) * (ATMTR - 0.07);
        }
        else
        {
            FRDIF = 1.0;
        }

        /* Diffuse photosynthetically active radiation parameter. */
        real DIFPP = FRDIF * ATMTR * 0.5 * SC;

        /* Populate and return the results structure. */
        result.DAYL = DAYL;
        result.DAYLP = DAYLP;
        result.SINLD = SINLD;
        result.COSLD = COSLD;
        result.DIFPP = DIFPP;
        result.ATMTR = ATMTR;
        result.DSINBE = DSINBE;
        result.ANGOT = ANGOT;

        return result;
    }

    /* ======================================================================= *
     * Penman (1948) Evaporation Method
     * ======================================================================= */

    /**
     * @brief Compute open water (E0), bare soil (ES0), and reference crop (ET0)
     *        evapotranspiration using the Penman (1948) combination equation.
     *
     * The Penman equation combines the energy balance (radiation term) with
     * the aerodynamic (wind and humidity) term:
     *
     *   ET = (delta * Rn + gamma * Ea) / (delta + gamma)
     *
     * where:
     *   delta  = slope of the saturation vapour pressure curve [kPa/deg C]
     *   Rn     = net radiation expressed as evaporative equivalent [mm/day]
     *   gamma  = psychrometric constant [kPa/deg C]
     *   Ea     = aerodynamic evaporative demand [mm/day]
     *
     * Three variants are computed with different albedo values:
     *   E0  (open water):     albedo = 0.05
     *   ES0 (bare soil):       albedo = 0.15
     *   ET0 (reference crop):  albedo = 0.25
     *
     * @param[in]  day     Date as std::tm (used to compute day-of-year)
     * @param[in]  LAT     Latitude [decimal degrees]
     * @param[in]  ELEV    Elevation [m a.s.l.]
     * @param[in]  TMIN    Daily minimum temperature [deg C]
     * @param[in]  TMAX    Daily maximum temperature [deg C]
     * @param[in]  AVRAD   Daily global solar radiation [J/m^2/d]
     * @param[in]  VAP     Daily vapour pressure [hPa]
     * @param[in]  WIND2   Wind speed at 2 m height [m/s]
     * @param[in]  ANGSTA  Angstrom A coefficient [dimensionless]
     * @param[in]  ANGSTB  Angstrom B coefficient [dimensionless]
     * @param[out] E0      Open water evaporation [mm/day]
     * @param[out] ES0     Bare soil evaporation [mm/day]
     * @param[out] ET0     Reference crop evapotranspiration [mm/day]
     */
    void penman(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                real AVRAD, real VAP, real WIND2, real ANGSTA, real ANGSTB,
                real &E0, real &ES0, real &ET0) const
    {
        /* --- Physical constants --- */
        real PSYCON = 0.67;              /* Psychrometric constant proportionality [mbar/deg C] */
        real REFCFW = 0.05;              /* Reflection coefficient (albedo) for open water [dimensionless] */
        real REFCFS = 0.15;              /* Reflection coefficient (albedo) for bare soil [dimensionless] */
        real REFCFC = 0.25;              /* Reflection coefficient (albedo) for reference crop [dimensionless] */
        real LHVAP = 2.45E6;             /* Latent heat of vaporization of water [J/kg]
                                          * (= 2.45 MJ/kg at ~20 deg C) */
        real STBC = 5.670373E-8 * 24 * 60 * 60;  /* Stefan-Boltzmann constant converted to [J/m^2/d/K^4]
                                                    * = 4.9E-3 approximately */

        /* --- Preparatory calculations --- */
        real TMPA = (TMIN + TMAX) / 2.0;  /* Mean daily temperature [deg C] */
        real TDIF = TMAX - TMIN;           /* Daily temperature range [deg C] */

        /* Wind function coefficient (Bu): depends on temperature range.
         * Larger diurnal range -> more turbulent mixing -> higher wind function. */
        real BU = 0.54 + 0.35 * limit(0.0, 1.0, (TDIF - 12.0) / 4.0);

        /* Barometric pressure [mbar], adjusted for altitude via the hypsometric equation. */
        real PBAR = 1013.0 * exp(-0.034 * ELEV / (TMPA + 273.0));

        /* Psychrometric constant [mbar/deg C]. */
        real GAMMA = PSYCON * PBAR / 1013.0;

        /* Saturated vapour pressure at mean temperature [mbar]. */
        real SVAP = 6.10588 * exp(17.32491 * TMPA / (TMPA + 238.102));

        /* Slope of the saturation vapour pressure curve at mean temperature [mbar/deg C]. */
        real DELTA = 238.102 * 17.32491 * SVAP / pow(TMPA + 238.102, 2);

        /* Ensure measured vapour pressure does not exceed saturation. */
        VAP = std::min(VAP, SVAP);

        /* Day of year for astronomical calculations. */
        int day_of_year = day.tm_yday + 1; /* tm_yday is 0-based */

        /* Atmospheric transmission and relative sunshine duration estimate. */
        AstroResults r = astro(day_of_year, LAT, AVRAD);
        real RELSSD = limit(0.0, 1.0, (r.ATMTR - abs(ANGSTA)) / abs(ANGSTB));

        /* --- Net outgoing long-wave radiation [J/m^2/d] ---
         * Depends on temperature, humidity, and cloudiness (via relative sunshine duration). */
        real RB = STBC * pow(TMPA + 273.0, 4) * (0.56 - 0.079 * sqrt(VAP)) * (0.1 + 0.9 * RELSSD);

        /* --- Net absorbed radiation converted to evaporative equivalent [mm/day] ---
         * Three variants based on different surface albedos. */
        real RNW = (AVRAD * (1.0 - REFCFW) - RB) / LHVAP;  /* Open water */
        real RNS = (AVRAD * (1.0 - REFCFS) - RB) / LHVAP;  /* Bare soil */
        real RNC = (AVRAD * (1.0 - REFCFC) - RB) / LHVAP;  /* Reference crop */

        /* --- Aerodynamic evaporative demand (drying power) [mm/day] ---
         * Two variants: one for open surfaces (E0, ES0) and one for crops (ET0)
         * with a higher wind function coefficient. */
        real EA = 0.26 * std::max(0.0, (SVAP - VAP)) * (0.5 + BU * WIND2);   /* Open/soil */
        real EAC = 0.26 * std::max(0.0, (SVAP - VAP)) * (1.0 + BU * WIND2);  /* Crop */

        /* --- Penman combination equation ---
         * E0  = open water evaporation
         * ES0 = bare soil evaporation
         * ET0 = reference crop evapotranspiration
         * All in [mm/day] */
        E0 = (DELTA * RNW + GAMMA * EA) / (DELTA + GAMMA);
        ES0 = (DELTA * RNS + GAMMA * EA) / (DELTA + GAMMA);
        ET0 = (DELTA * RNC + GAMMA * EAC) / (DELTA + GAMMA);

        /* Ensure non-negative reference evaporation. */
        E0 = std::max(0.0, E0);
        ES0 = std::max(0.0, ES0);
        ET0 = std::max(0.0, ET0);
    }

    /* ======================================================================= *
     * Penman-Monteith Reference Evapotranspiration Method
     * ======================================================================= */

    /**
     * @brief Compute reference crop evapotranspiration using the Penman-Monteith
     *        equation (FAO-56 style formulation).
     *
     * This is the physically-based combination equation that accounts for
     * both energy availability and aerodynamic transport, with explicit
     * surface resistance:
     *
     *   ET0 = [delta * (Rn - G) + gamma * (900/(T+273)) * u2 * (es - ea)]
     *         / [delta + gamma * (1 + rs/ra)]
     *
     * where:
     *   delta = slope of saturation vapour pressure curve [kPa/deg C]
     *   Rn    = net radiation [mm/day]
     *   G     = soil heat flux [mm/day] (assumed 0 for daily time step)
     *   gamma = psychrometric constant [kPa/deg C]
     *   u2    = wind speed at 2 m height [m/s]
     *   es    = mean saturation vapour pressure [kPa]
     *   ea    = actual vapour pressure [kPa]
     *   rs    = surface resistance [s/m] (70 s/m for reference crop)
     *   ra    = aerodynamic resistance [s/m]
     *
     * @param[in]  day     Date as std::tm (used for day-of-year calculation)
     * @param[in]  LAT     Site latitude [decimal degrees]
     * @param[in]  ELEV    Site elevation [m a.s.l.]
     * @param[in]  TMIN    Daily minimum temperature [deg C]
     * @param[in]  TMAX    Daily maximum temperature [deg C]
     * @param[in]  AVRAD   Daily global solar radiation [J/m^2/d]
     * @param[in]  VAP     Daily vapour pressure [hPa]
     * @param[in]  WIND2   Wind speed at 2 m height [m/s]
     * @return Reference crop evapotranspiration [mm/day]
     *
     * Typical ET0 ranges (temperate climate):
     *   Winter: 0.5 - 2.0 mm/day
     *   Summer: 3.0 - 7.0 mm/day
     */
    real penman_monteith(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                         real AVRAD, real VAP, real WIND2) const
    {
        /* --- Physical constants --- */
        real PSYCON = 0.665;             /* Psychrometric constant proportionality [kPa/deg C] */
        real REFCFC = 0.23;              /* Albedo for reference crop (grass) [dimensionless] */
        real CRES = 70.0;                /* Surface (canopy) resistance [s/m]
                                          * 70 s/m is the FAO-56 standard for a well-watered
                                          * reference grass crop of height 0.12 m */
        real LHVAP = 2.45E6;             /* Latent heat of vaporization [J/kg] */
        real STBC = 4.903E-3;            /* Stefan-Boltzmann constant [J/m^2/d/K^4]
                                          * (5.670373E-8 W/m^2/K^4 * 86400 s/d) */
        real G = 0.0;                    /* Soil heat flux density [mm/day]
                                          * Assumed zero for daily time steps (FAO-56) */

        /* Mean daily temperature [deg C]. */
        real TMPA = (TMIN + TMAX) / 2.0;

        /* Convert vapour pressure from hPa to kPa. */
        real VAP_kPa = hPa2kPa(VAP);

        /* Atmospheric pressure at the site elevation [kPa].
         * Uses the barometric formula with a standard lapse rate of 0.0065 K/m
         * and a reference temperature of 293 K (20 deg C). */
        real T = 293.0;
        real PATM = 101.3 * pow((T - (0.0065 * ELEV)) / T, 5.26);

        /* Psychrometric constant [kPa/deg C]. */
        real GAMMA = PSYCON * PATM * 1.0E-3;

        /* Saturated vapour pressure at mean temperature [kPa]. */
        real SVAP_TMPA = SatVapourPressure(TMPA);

        /* Slope of the saturation vapour pressure curve at mean temperature [kPa/deg C]. */
        real DELTA = (4098.0 * SVAP_TMPA) / pow(TMPA + 237.3, 2);

        /* Mean daily saturated vapour pressure [kPa], computed as the average
         * of saturation pressures at Tmin and Tmax (FAO-56 recommendation). */
        real SVAP_TMAX = SatVapourPressure(TMAX);
        real SVAP_TMIN = SatVapourPressure(TMIN);
        real SVAP = (SVAP_TMAX + SVAP_TMIN) / 2.0;

        /* Ensure measured vapour pressure does not exceed saturation. */
        VAP_kPa = std::min(VAP_kPa, SVAP);

        /* --- Net outgoing long-wave radiation ---
         * Computed in two steps: temperature-dependent emission, then
         * cloudiness correction. */
        /* Step 1: Black-body emission at Tmin, Tmax, averaged [J/m^2/d]. */
        real STB_TMAX = STBC * pow(Celsius2Kelvin(TMAX), 4);
        real STB_TMIN = STBC * pow(Celsius2Kelvin(TMIN), 4);
        real RNL_TMP = ((STB_TMAX + STB_TMIN) / 2.0) * (0.34 - 0.14 * sqrt(VAP_kPa));

        /* Step 2: Clear-sky radiation from Angot TOA radiation [J/m^2/d].
         * The clear-sky transmissivity increases with elevation due to
         * reduced atmospheric path length. */
        int day_of_year = day.tm_yday + 1;
        AstroResults r = astro(day_of_year, LAT, AVRAD);
        real CSKYRAD = (0.75 + (2e-05 * ELEV)) * r.ANGOT;

        real ET0 = 0.0;
        if (CSKYRAD > 0)
        {
            /* Cloudiness correction for long-wave radiation. */
            real RNL = RNL_TMP * (1.35 * (AVRAD / CSKYRAD) - 0.35);

            /* Net radiation converted to evaporative equivalent [mm/day]. */
            real RN = ((1.0 - REFCFC) * AVRAD - RNL) / LHVAP;

            /* Aerodynamic evaporation term [mm/day].
             * Combines the vapour pressure deficit with wind speed and
             * a temperature-dependent scaling factor (900/(T+273)). */
            real EA = ((900.0 / (TMPA + 273.0)) * WIND2 * (SVAP - VAP_kPa));

            /* Modified psychrometric constant [kPa/deg C].
             * Incorporates the surface and aerodynamic resistances:
             *   gamma* = gamma * (1 + rs / (208 * u2))
             * where rs = CRES (70 s/m) and 208 = ra_multiplier. */
            real MGAMMA = GAMMA * (1.0 + (CRES / 208.0 * WIND2));

            /* Penman-Monteith combination equation [mm/day]. */
            ET0 = (DELTA * (RN - G)) / (DELTA + MGAMMA) + (GAMMA * EA) / (DELTA + MGAMMA);

            /* Ensure non-negative result. */
            ET0 = std::max(0.0, ET0);
        }
        return ET0;
    }

    /* ======================================================================= *
     * Combined ET Calculation Entry Point
     * ======================================================================= */

    /**
     * @brief Compute all three reference evapotranspiration rates for a single day.
     *
     * This method delegates to:
     *   - penman() for E0 (open water) and ES0 (bare soil) using the Penman (1948) method
     *   - penman_monteith() for ET0 (reference crop) using the FAO-56 Penman-Monteith method
     *
     * @param[in]  day     Date as std::tm
     * @param[in]  LAT     Site latitude [decimal degrees]
     * @param[in]  ELEV    Site elevation [m a.s.l.]
     * @param[in]  TMIN    Daily minimum temperature [deg C]
     * @param[in]  TMAX    Daily maximum temperature [deg C]
     * @param[in]  AVRAD   Daily global solar radiation [J/m^2/d]
     * @param[in]  VAP     Daily vapour pressure [hPa]
     * @param[in]  WIND    Wind speed at 2 m height [m/s]
     * @param[in]  ANGSTA  Angstrom A coefficient [dimensionless]
     * @param[in]  ANGSTB  Angstrom B coefficient [dimensionless]
     * @param[out] E0      Open water evaporation [mm/day]
     * @param[out] ES0     Bare soil evaporation [mm/day]
     * @param[out] ET0     Reference crop evapotranspiration [mm/day]
     */
    void calculateET(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                     real AVRAD, real VAP, real WIND, real ANGSTA, real ANGSTB,
                     real &E0, real &ES0, real &ET0) const
    {
        /* Compute E0 and ES0 using the Penman (1948) method.
         * Also computes a Penman-style ET0 (et0_penman) which is discarded
         * in favour of the Penman-Monteith result below. */
        real et0_penman;
        penman(day, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND, ANGSTA, ANGSTB,
               E0, ES0, et0_penman);

        /* Compute ET0 using the physically-based Penman-Monteith method,
         * which accounts for surface and aerodynamic resistances explicitly. */
        ET0 = penman_monteith(day, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND);
    }
};

#endif
