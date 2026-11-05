#ifndef _SERGHEI_MPI_UTILS_H_
#define _SERGHEI_MPI_UTILS_H_

#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>

// All model components currently share MPI_COMM_WORLD.  Centralising error
// handling prevents a failed nonblocking operation from being overwritten by
// the return code of a later MPI call.
inline void sergheiMpiCheck(const int code, const char *operation,
                            const char *file, const int line,
                            MPI_Comm comm = MPI_COMM_WORLD)
{
    if (code == MPI_SUCCESS)
        return;

    int rank = -1;
    int initialized = 0;
    int finalized = 0;
    MPI_Initialized(&initialized);
    if (initialized)
        MPI_Finalized(&finalized);
    MPI_Comm error_comm = (comm == MPI_COMM_NULL) ? MPI_COMM_WORLD : comm;
    char message[MPI_MAX_ERROR_STRING] = "MPI failure before initialization";
    int length = 33;
    if (initialized && !finalized)
    {
        MPI_Error_string(code, message, &length);
        MPI_Comm_rank(error_comm, &rank);
    }
    std::cerr << "[MPI ERROR] rank=" << rank << " operation=" << operation
              << " at " << file << ":" << line << ": "
              << std::string(message, length) << std::endl;
    if (initialized && !finalized)
        MPI_Abort(error_comm, code);
    std::abort();
}

#define SERGHEI_MPI_CHECK_COMM(call, comm) \
    sergheiMpiCheck((call), #call, __FILE__, __LINE__, (comm))

#define SERGHEI_MPI_CHECK(call) \
    SERGHEI_MPI_CHECK_COMM((call), MPI_COMM_WORLD)

inline unsigned long long sergheiFileFingerprint(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    unsigned long long hash = 1469598103934665603ULL;
    if (!input.is_open())
        return hash ^ 0x9e3779b97f4a7c15ULL;
    char ch;
    while (input.get(ch))
    {
        hash ^= static_cast<unsigned char>(ch);
        hash *= 1099511628211ULL;
    }
    return hash;
}

inline bool sergheiVerifyFileConsistency(const std::string &path, int rank)
{
    const unsigned long long local = sergheiFileFingerprint(path);
    unsigned long long min_hash = 0, max_hash = 0;
    SERGHEI_MPI_CHECK(MPI_Allreduce(&local, &min_hash, 1, MPI_UNSIGNED_LONG_LONG,
                                    MPI_MIN, MPI_COMM_WORLD));
    SERGHEI_MPI_CHECK(MPI_Allreduce(&local, &max_hash, 1, MPI_UNSIGNED_LONG_LONG,
                                    MPI_MAX, MPI_COMM_WORLD));
    if (min_hash == max_hash)
        return true;
    std::cerr << "[RT MPI ERROR] input differs across ranks: " << path
              << " (rank " << rank << ", fingerprint=" << local << ")"
              << std::endl;
    return false;
}

#endif
