Remotes globals_to_linear_owners(Read<GO> globals, GO total, I32 comm_size);
Remotes globals_to_linear_owners(CommPtr comm, Read<GO> globals, GO total);
LO linear_partition_size(GO total, I32 comm_size, I32 comm_rank);
LO linear_partition_size(CommPtr comm, GO total);
GO find_total_globals(CommPtr comm, Read<GO> globals);
