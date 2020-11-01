#pragma once
#include <mpi.h>

class MpiHolder
{
public:
    MpiHolder(int* argc, char*** argv)
    {
        MPI_Init(argc, argv);
    }
    ~MpiHolder()
    {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
};