/* see hilbert.cpp */
namespace hilbert {
/* original code had unsigned int here, we are going to 64 bits */
typedef std::uint64_t coord_t; //char,short,int for up to 8,16,32 bits per word
void TransposetoAxes( coord_t* X, int b, int n ); // position, #bits, dimension
void AxestoTranspose( coord_t* X, int b, int n ); // position, #bits, dimension
}
