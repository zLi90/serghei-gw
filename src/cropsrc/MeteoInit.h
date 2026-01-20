/* -*- mode: c++; c-default-style: "linux" -*- */

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

class MeteoInit
{

public:
    // Main function to read meteo data
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
        // Temporary storage on Host using std::vector
        std::vector<std::vector<real>> raw_data;

        // Default location parameters (Placeholder logic, adapt as needed)
        ms.lat = -999.0;
        ms.lon = -999.0;
        ms.elev = -999.0;
        // Default Angstrom coefficients (typical values for clear sky conditions)
        ms.angstA = 0.18;
        ms.angstB = 0.50;

        // Store the file content for parsing location info
        std::vector<std::string> file_lines;

        while (std::getline(fInStream, line))
        {
            file_lines.push_back(line);

            // Trim whitespace
            size_t first = line.find_first_not_of(" \t\r\n");
            if (first == std::string::npos)
                continue; // Empty line

            // Check for location info first (before data parsing)
            parseLocationInfo(line, ms);

            // Skip comments for data parsing
            if (line[first] == '/' && line.size() > first + 1 && line[first + 1] == '/')
            {
                continue;
            }

            // Data lines (start with digit or minus sign)
            if (isdigit(line[first]) || line[first] == '-')
            {
                std::stringstream ss(line);
                std::vector<real> row_vals;
                real val;
                while (ss >> val)
                {
                    row_vals.push_back(val);
                }

                // Expected columns: Time, Year, Day, Irrad, Tmin, Tmax, Vap, Wind, Rain
                // Total 9 columns minimum
                if (row_vals.size() >= 9)
                {
                    raw_data.push_back(row_vals);
                }
            }
        }
        fInStream.close();

        size_t n = raw_data.size();
        if (n == 0)
        {
            if (par.masterproc)
                std::cerr << RERROR "No valid meteo records found.\n";
            return 0;
        }

        // ================== NEW CODE START ==================
        // 1. Extract Start Date from the first record
        // Assuming Column 1 is Year, Column 2 is Day (DOY)
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
            // Fallback default
            ms.start_year = 2000;
            ms.start_doy = 1;
        }
        // ================== NEW CODE END ==================

        // 1. Allocate Kokkos Views on Device
        ms.allocate(n);

        // 2. Create Host Mirrors (CPU accessible views)
        auto h_time = Kokkos::create_mirror_view(ms.time);
        auto h_irrad = Kokkos::create_mirror_view(ms.irrad);
        auto h_tmin = Kokkos::create_mirror_view(ms.tmin);
        auto h_tmax = Kokkos::create_mirror_view(ms.tmax);
        auto h_vap = Kokkos::create_mirror_view(ms.vap);
        auto h_wind = Kokkos::create_mirror_view(ms.wind);
        auto h_rain = Kokkos::create_mirror_view(ms.rain);

        // Mirrors for derived ET variables (initialized to 0)
        auto h_e0 = Kokkos::create_mirror_view(ms.e0);
        auto h_es0 = Kokkos::create_mirror_view(ms.es0);
        auto h_et0 = Kokkos::create_mirror_view(ms.et0);

        // 3. Fill Host Mirrors with data from vector
        for (size_t i = 0; i < n; ++i)
        {
            // Mapping based on: Time, Year, Day, Irrad, Tmin, Tmax, Vap, Wind, Rain
            // 0     1     2    3      4     5     6    7     8
            h_time(i) = raw_data[i][0];
            h_irrad(i) = raw_data[i][3];
            h_tmin(i) = raw_data[i][4];
            h_tmax(i) = raw_data[i][5];
            h_vap(i) = raw_data[i][6];
            h_wind(i) = raw_data[i][7];
            h_rain(i) = raw_data[i][8];

            // Calculate ET values using proper Penman and Penman-Monteith methods
            real TMIN = raw_data[i][4];
            real TMAX = raw_data[i][5];
            real IRRAD = raw_data[i][3]; // J/m2/d
            real VAP = raw_data[i][6];   // hPa
            real WIND = raw_data[i][7];  // m/s

            // ================== MODIFIED DATE LOGIC START ==================
            // Use the actual Year and DOY from the file for ET calculation
            int current_year = (int)raw_data[i][1];
            int current_doy = (int)raw_data[i][2];

            std::tm date_tm = {};
            date_tm.tm_year = current_year - 1900;
            date_tm.tm_mday = 1;
            date_tm.tm_mon = 0;    // Jan 1st
            std::mktime(&date_tm); // Normalize
            // Add DOY offset (tm_yday is 0-based, so adding DOY-1 days)
            date_tm.tm_mday += (current_doy - 1);
            std::mktime(&date_tm); // Re-normalize to get correct month/day
            // ================== MODIFIED DATE LOGIC END ==================

            // Calculate E0, ES0, ET0 using Penman and Penman-Monteith
            real e0, es0, et0;
            calculateET(date_tm, ms.lat, ms.elev, TMIN, TMAX, IRRAD, VAP, WIND,
                        ms.angstA, ms.angstB, e0, es0, et0);

            // Convert from mm/day to cm/day (1 mm = 0.1 cm)
            // h_e0(i)  = e0 * 0.1;   // cm/day
            // h_es0(i) = es0 * 0.1;  // cm/day
            // ms.et0(i) = et0 * 0.1;  // cm/day

            ms.e0(i) = e0 * 0.1;   // cm/day
            ms.es0(i) = es0 * 0.1; // cm/day
            ms.et0(i) = et0 * 0.1; // cm/day
#if DEBUG_CROP_GROWTH_MODEL
            printf("MeteoInit: i=%zu, DOY=%d, E0=%.6f cm/day, ES0=%.6f cm/day, ET0=%.6f cm/day\n",
                   i, day_of_year, h_e0(i), h_es0(i), ms.et0(i));
#endif
        }

        // 4. Deep Copy from Host Mirrors to Device Views
        Kokkos::deep_copy(ms.time, h_time);
        Kokkos::deep_copy(ms.irrad, h_irrad);
        Kokkos::deep_copy(ms.tmin, h_tmin);
        Kokkos::deep_copy(ms.tmax, h_tmax);
        Kokkos::deep_copy(ms.vap, h_vap);
        Kokkos::deep_copy(ms.wind, h_wind);
        Kokkos::deep_copy(ms.rain, h_rain);

        // Kokkos::deep_copy(ms.e0, h_e0);
        // Kokkos::deep_copy(ms.es0, h_es0);
        // Kokkos::deep_copy(ms.et0, h_et0);

#if DEBUG_CROP_GROWTH_MODEL
        //  Copy device view back to host and print all ET0 values
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
    // Helper function to parse location information from lines
    void parseLocationInfo(const std::string &line, MeteoState &ms)
    {
        std::string trimmed = line;

        // Remove leading // and whitespace if present
        size_t pos = trimmed.find("//");
        if (pos != std::string::npos)
        {
            trimmed = trimmed.substr(pos + 2);
        }

        // Trim leading/trailing whitespace
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

        // Parse Lat: value
        if (trimmed.find("Lat:") != std::string::npos)
        {
            size_t lat_pos = trimmed.find("Lat:");
            std::string substr = trimmed.substr(lat_pos + 4);
            // Remove spaces and commas
            substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
            substr.erase(std::remove(substr.begin(), substr.end(), ','), substr.end());
            // Extract number
            std::stringstream ss(substr);
            double val;
            if (ss >> val)
            {
                ms.lat = val;
            }
        }

        // Parse Lon: value
        if (trimmed.find("Lon:") != std::string::npos)
        {
            size_t lon_pos = trimmed.find("Lon:");
            std::string substr = trimmed.substr(lon_pos + 4);
            // Remove spaces and commas
            substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
            substr.erase(std::remove(substr.begin(), substr.end(), ','), substr.end());
            // Extract number
            std::stringstream ss(substr);
            double val;
            if (ss >> val)
            {
                ms.lon = val;
            }
        }

        // Parse Elev: value
        if (trimmed.find("Elev:") != std::string::npos)
        {
            size_t elev_pos = trimmed.find("Elev:");
            std::string substr = trimmed.substr(elev_pos + 5);
            // Remove spaces and commas
            substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
            substr.erase(std::remove(substr.begin(), substr.end(), ','), substr.end());
            // Extract number
            std::stringstream ss(substr);
            double val;
            if (ss >> val)
            {
                ms.elev = val;
            }
        }

        // Parse Angstrom A and B coefficients
        if (trimmed.find("AngstromA:") != std::string::npos)
        {
            size_t anga_pos = trimmed.find("AngstromA:");
            std::string substr = trimmed.substr(anga_pos + 11);
            substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
            substr.erase(std::remove(substr.begin(), substr.end(), ','), substr.end());
            std::stringstream ss(substr);
            double val;
            if (ss >> val)
            {
                ms.angstA = val;
            }
        }

        if (trimmed.find("AngstromB:") != std::string::npos)
        {
            size_t angb_pos = trimmed.find("AngstromB:");
            std::string substr = trimmed.substr(angb_pos + 11);
            substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
            substr.erase(std::remove(substr.begin(), substr.end(), ','), substr.end());
            std::stringstream ss(substr);
            double val;
            if (ss >> val)
            {
                ms.angstB = val;
            }
        }
    }

    // Helper functions for ET calculations
    real hPa2kPa(real hPa) const
    {
        return hPa / 10.0;
    }

    real Celsius2Kelvin(real celsius) const
    {
        return celsius + 273.16;
    }

    real SatVapourPressure(real temp) const
    {
        return 0.6108 * exp((17.27 * temp) / (237.3 + temp));
    }

    real limit(real vmin, real vmax, real v) const
    {
        if (v < vmin)
            return vmin;
        if (v > vmax)
            return vmax;
        return v;
    }

    // Astro calculation for day length and radiation
    struct AstroResults
    {
        real DAYL, DAYLP, SINLD, COSLD, DIFPP, ATMTR, DSINBE, ANGOT;
    };

    AstroResults astro(int day_of_year, real latitude, real radiation) const
    {
        AstroResults result = {};

        if (abs(latitude) > 90.0)
        {
            return result;
        }

        real RAD = M_PI / 180.0;
        real LAT = latitude;

        // Declination
        real DEC = -asin(sin(23.45 * RAD) * cos(2.0 * M_PI * (day_of_year + 10.0) / 365.0));
        real SC = 1370.0 * (1.0 + 0.033 * cos(2.0 * M_PI * day_of_year / 365.0));

        // Intermediate variables
        real SINLD = sin(RAD * LAT) * sin(DEC);
        real COSLD = cos(RAD * LAT) * cos(DEC);
        real AOB = SINLD / COSLD;

        // Day length for base = 0 degrees
        real DAYL = 0.0;
        real DSINB = 0.0;
        real DSINBE = 0.0;

        if (abs(AOB) <= 1.0)
        {
            DAYL = 12.0 * (1.0 + 2.0 * asin(AOB) / M_PI);
            DSINB = 3600.0 * (DAYL * SINLD + 24.0 * COSLD * sqrt(1.0 - AOB * AOB) / M_PI);
            DSINBE = 3600.0 * (DAYL * (SINLD + 0.4 * (SINLD * SINLD + COSLD * COSLD * 0.5)) +
                               12.0 * COSLD * (2.0 + 3.0 * 0.4 * SINLD) * sqrt(1.0 - AOB * AOB) / M_PI);
        }
        else
        {
            if (AOB > 1.0)
                DAYL = 24.0;
            if (AOB < -1.0)
                DAYL = 0.0;
            DSINB = 3600.0 * (DAYL * SINLD);
            DSINBE = 3600.0 * (DAYL * (SINLD + 0.4 * (SINLD * SINLD + COSLD * COSLD * 0.5)));
        }

        // Day length for base = -4 degrees
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

        // Extraterrestrial radiation and atmospheric transmission
        real ANGOT = SC * DSINB;
        real ATMTR = 0.0;
        if (DAYL > 0.0)
        {
            ATMTR = radiation / ANGOT;
        }

        // Estimate fraction diffuse irradiation
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

        real DIFPP = FRDIF * ATMTR * 0.5 * SC;

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

    // Penman method for E0, ES0, ET0
    void penman(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                real AVRAD, real VAP, real WIND2, real ANGSTA, real ANGSTB,
                real &E0, real &ES0, real &ET0) const
    {

        // Constants
        real PSYCON = 0.67;
        real REFCFW = 0.05;
        real REFCFS = 0.15;
        real REFCFC = 0.25;
        real LHVAP = 2.45E6;
        real STBC = 5.670373E-8 * 24 * 60 * 60; // =4.9E-3

        // Preparatory calculations
        real TMPA = (TMIN + TMAX) / 2.0;
        real TDIF = TMAX - TMIN;
        real BU = 0.54 + 0.35 * limit(0.0, 1.0, (TDIF - 12.0) / 4.0);

        // Barometric pressure (mbar) and psychrometric constant (mbar/Celsius)
        real PBAR = 1013.0 * exp(-0.034 * ELEV / (TMPA + 273.0));
        real GAMMA = PSYCON * PBAR / 1013.0;

        // Saturated vapour pressure and its derivative
        real SVAP = 6.10588 * exp(17.32491 * TMPA / (TMPA + 238.102));
        real DELTA = 238.102 * 17.32491 * SVAP / pow(TMPA + 238.102, 2);
        VAP = std::min(VAP, SVAP);

        // Get day of year
        int day_of_year = day.tm_yday + 1; // tm_yday is 0-based

        // Calculate atmospheric transmission
        AstroResults r = astro(day_of_year, LAT, AVRAD);
        real RELSSD = limit(0.0, 1.0, (r.ATMTR - abs(ANGSTA)) / abs(ANGSTB));

        // Net outgoing long-wave radiation (J/m2/d)
        real RB = STBC * pow(TMPA + 273.0, 4) * (0.56 - 0.079 * sqrt(VAP)) * (0.1 + 0.9 * RELSSD);

        // Net absorbed radiation, expressed in mm/d
        real RNW = (AVRAD * (1.0 - REFCFW) - RB) / LHVAP;
        real RNS = (AVRAD * (1.0 - REFCFS) - RB) / LHVAP;
        real RNC = (AVRAD * (1.0 - REFCFC) - RB) / LHVAP;

        // Evaporative demand of the atmosphere (mm/d)
        real EA = 0.26 * std::max(0.0, (SVAP - VAP)) * (0.5 + BU * WIND2);
        real EAC = 0.26 * std::max(0.0, (SVAP - VAP)) * (1.0 + BU * WIND2);

        // Penman formula (1948)
        E0 = (DELTA * RNW + GAMMA * EA) / (DELTA + GAMMA);
        ES0 = (DELTA * RNS + GAMMA * EA) / (DELTA + GAMMA);
        ET0 = (DELTA * RNC + GAMMA * EAC) / (DELTA + GAMMA);

        // Ensure reference evaporation >= 0
        E0 = std::max(0.0, E0);
        ES0 = std::max(0.0, ES0);
        ET0 = std::max(0.0, ET0);
    }

    // Penman-Monteith method for ET0
    real penman_monteith(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                         real AVRAD, real VAP, real WIND2) const
    {

        // Constants
        real PSYCON = 0.665;
        real REFCFC = 0.23;
        real CRES = 70.0;
        real LHVAP = 2.45E6;
        real STBC = 4.903E-3;
        real G = 0.0; // Soil heat flux

        // Mean daily temperature
        real TMPA = (TMIN + TMAX) / 2.0;

        // Vapour pressure to kPa
        real VAP_kPa = hPa2kPa(VAP);

        // Atmospheric pressure at standard temperature of 293K (kPa)
        real T = 293.0;
        real PATM = 101.3 * pow((T - (0.0065 * ELEV)) / T, 5.26);

        // Psychrometric constant (kPa/Celsius)
        real GAMMA = PSYCON * PATM * 1.0E-3;

        // Derivative of SVAP with respect to mean temperature
        real SVAP_TMPA = SatVapourPressure(TMPA);
        real DELTA = (4098.0 * SVAP_TMPA) / pow(TMPA + 237.3, 2);

        // Daily average saturated vapour pressure [kPa] from min/max temperature
        real SVAP_TMAX = SatVapourPressure(TMAX);
        real SVAP_TMIN = SatVapourPressure(TMIN);
        real SVAP = (SVAP_TMAX + SVAP_TMIN) / 2.0;

        // Measured vapour pressure not to exceed saturated vapour pressure
        VAP_kPa = std::min(VAP_kPa, SVAP);

        // Longwave radiation at Tmax, Tmin (J/m2/d)
        real STB_TMAX = STBC * pow(Celsius2Kelvin(TMAX), 4);
        real STB_TMIN = STBC * pow(Celsius2Kelvin(TMIN), 4);
        real RNL_TMP = ((STB_TMAX + STB_TMIN) / 2.0) * (0.34 - 0.14 * sqrt(VAP_kPa));

        // Clear Sky radiation [J/m2/DAY] from Angot TOA radiation
        int day_of_year = day.tm_yday + 1;
        AstroResults r = astro(day_of_year, LAT, AVRAD);
        real CSKYRAD = (0.75 + (2e-05 * ELEV)) * r.ANGOT;

        real ET0 = 0.0;
        if (CSKYRAD > 0)
        {
            // Final net outgoing longwave radiation [J/m2/day]
            real RNL = RNL_TMP * (1.35 * (AVRAD / CSKYRAD) - 0.35);

            // Radiative evaporation equivalent for the reference surface [mm/DAY]
            real RN = ((1.0 - REFCFC) * AVRAD - RNL) / LHVAP;

            // Aerodynamic evaporation equivalent [mm/day]
            real EA = ((900.0 / (TMPA + 273.0)) * WIND2 * (SVAP - VAP_kPa));

            // Modified psychometric constant (gamma*)[kPa/C]
            real MGAMMA = GAMMA * (1.0 + (CRES / 208.0 * WIND2));

            // Reference ET in mm/day
            ET0 = (DELTA * (RN - G)) / (DELTA + MGAMMA) + (GAMMA * EA) / (DELTA + MGAMMA);

            ET0 = std::max(0.0, ET0);
        }
        return ET0;
    }

    // Main ET calculation function
    void calculateET(const std::tm &day, real LAT, real ELEV, real TMIN, real TMAX,
                     real AVRAD, real VAP, real WIND, real ANGSTA, real ANGSTB,
                     real &E0, real &ES0, real &ET0) const
    {

        // Calculate using Penman method for E0 and ES0
        real et0_penman;
        penman(day, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND, ANGSTA, ANGSTB,
               E0, ES0, et0_penman);

        // Calculate ET0 using Penman-Monteith method
        ET0 = penman_monteith(day, LAT, ELEV, TMIN, TMAX, AVRAD, VAP, WIND);
    }
};

#endif
