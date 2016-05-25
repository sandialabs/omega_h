/* A reproducible sum of floating-point values.
   this operation is one of the key places where
   a program's output begins to depend on parallel
   partitioning and traversal order, because
   floating-point values do not produce the same
   sum when added in a different order.

   IEEE 754 64-bit floating point format is assumed,
   which has 52 bits in the fraction.

   The idea here is to add the numbers as fixed-point values.
   max_exponent() finds the largest exponent (e) such that
   all values are (<= 2^(e)).
   We then use the value (2^(e - 52)) as the unit, and sum all
   values as integers in that unit.
   This is guaranteed to be at least as accurate as the
   worst-case ordering of the values, i.e. being added
   in order of decreasing magnitude.

   If we used a 64-bit integer type, we would only be
   able to reliably add up to (2^12 = 4096) values
   (64 - 52 = 12).
   Thus we use a 128-bit integer type.
   This allows us to reliably add up to (2^76 > 10^22)
   values. By comparison, supercomputers today
   support a maximum of one million MPI ranks (10^6)
   and each rank typically can't hold more than
   one billion values (10^9), for a total of (10^15) values.
*/

int max_exponent(Reals a);
Real repro_sum(Reals a);
Real repro_sum(CommPtr comm, Reals a);
