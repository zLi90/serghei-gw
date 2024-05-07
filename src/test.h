inline real calculate_cfl_dt(const realArr2 &velocity, real dx, real dy, real dz, real cfl)
{
    real max_velocity = 0.0;
    Kokkos::parallel_reduce(
        "MaxVelocity", gdom.nCellMem, KOKKOS_LAMBDA(int iGlob, real &max_vel) {
            real u = velocity(iGlob, 0);
            real v = velocity(iGlob, 1);
            real w = velocity(iGlob, 2);
            real vel = sqrt(u * u + v * v + w * w);
            max_vel = fmax(max_vel, vel);
        },
        Kokkos::Max<real>(max_velocity));

    // Ensure denominator is not zero
    max_velocity += 1e-10;

    // Calculate CFL timestep
    real dt_cfl = cfl * fmin(fmin(dx, dy), dz) / max_velocity;

    return dt_cfl;
}

void rt_solve(RTState &rt, GwState &gw, GwDomain &gdom, GwMPI &gmpi, Parallel &par)
{
    real cfl = 0.5; // CFL number
    real dx = gdom.dx;
    real dy = gdom.dy;
    real dz = gdom.dz;

    // Calculate CFL timestep
    real dt_cfl = calculate_cfl_dt(gw.velocity, dx, dy, dz, cfl);

    // Iterate over possible timestep values
    for (real dt = dt_cfl; dt > 1e-10; dt *= 0.5)
    {
        // Solve transport equation with current timestep
        advection_dispersion(rt, gw, gdom, dt);

        // Check for convergence
        if (check_convergence())
        {
            // If converged, update current timestep and exit loop
            rt.dt = dt;
            break;
        }
    }

    // Update concentration arrays
    Kokkos::parallel_for(
        gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            rt.c(iGlob, 0) = rt.c(iGlob, 1);
        });

    // MPI communication
    gmpi.mpi_sendrecv(rt.c, gdom, par);

    // Write output for reactive transport
    // io.outputIniRT(outFolder);
}

inline void advection_dispersion(RTState &rt, GwState &gw, GwDomain &gdom, real dt)
{
    real dx = gdom.dx;
    real dy = gdom.dy;
    real dz = gdom.dz;

    Kokkos::parallel_for(
        gdom.nCellMem, KOKKOS_LAMBDA(int iGlob) {
            // Calculate advection contribution
            real u = gw.velocity(iGlob, 0);
            real v = gw.velocity(iGlob, 1);
            real w = gw.velocity(iGlob, 2);
            real advection = -u * (rt.c(iGlob, 1) - rt.c(iGlob - gdom.nx * gdom.ny, 1)) / dx - v * (rt.c(iGlob, 1) - rt.c(iGlob - gdom.nx, 1)) / dy - w * (rt.c(iGlob, 1) - rt.c(iGlob - 1, 1)) / dz;

            // Calculate dispersion contribution
            // (To be implemented)

            // Update concentration
            rt.c(iGlob, 0) = rt.c(iGlob, 1) + dt * advection;
        });
}
