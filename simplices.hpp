#define DIMS 4

enum {
  VERT,
  EDGE,
  FACE,
  TRI = FACE,
  REGION,
  TET = REGION
};

CONSTANT static I8 const fev0[] = {0,1};
CONSTANT static I8 const fev1[] = {1,2};
CONSTANT static I8 const fev2[] = {2,0};
CONSTANT static I8 const* const fev_[] = {
  fev0,fev1,fev2};
CONSTANT static I8 const* const* const f_v_[] = {
  0,fev_};
CONSTANT static I8 const rev0[] = {0,1};
CONSTANT static I8 const rev1[] = {1,2};
CONSTANT static I8 const rev2[] = {2,0};
CONSTANT static I8 const rev3[] = {0,3};
CONSTANT static I8 const rev4[] = {1,3};
CONSTANT static I8 const rev5[] = {2,3};
CONSTANT static I8 const* const rev_[] = {
  rev0,rev1,rev2,rev3,rev4,rev5};
CONSTANT static I8 const rfv0[] = {0,2,1};
CONSTANT static I8 const rfv1[] = {0,1,3};
CONSTANT static I8 const rfv2[] = {1,2,3};
CONSTANT static I8 const rfv3[] = {2,0,3};
CONSTANT static I8 const* const rfv_[] = {
  rfv0,rfv1,rfv2,rfv3};
CONSTANT static I8 const* const* const r_v_[] = {
  0,rev_,rfv_};
CONSTANT static I8 const* const* const* const down_templates[] = {
  0,0,f_v_,r_v_};

struct UpTemplate { I8 up; I8 which_down; bool is_flipped; };
CONSTANT static UpTemplate const fve0[] = {{0,0,0},{2,1,0}};
CONSTANT static UpTemplate const fve1[] = {{1,0,0},{0,1,0}};
CONSTANT static UpTemplate const fve2[] = {{2,0,0},{1,1,0}};
CONSTANT static UpTemplate const* const fve_[] = {
  fve0,fve1,fve2};
CONSTANT static UpTemplate const* const* const f_u_[] = {
  fve_};
CONSTANT static UpTemplate const rve0[] = {{0,0,0},{2,1,0},{3,0,0}};
CONSTANT static UpTemplate const rve1[] = {{1,0,0},{0,1,0},{4,0,0}};
CONSTANT static UpTemplate const rve2[] = {{2,0,0},{1,1,0},{5,0,0}};
CONSTANT static UpTemplate const rve3[] = {{3,1,0},{4,1,0},{5,1,0}};
CONSTANT static UpTemplate const* const rve_[] = {
  rve0,rve1,rve2,rve3};
CONSTANT static UpTemplate const ref0[] = {{0,2,1},{1,0,0}};
CONSTANT static UpTemplate const ref1[] = {{0,1,1},{2,0,0}};
CONSTANT static UpTemplate const ref2[] = {{0,0,1},{3,0,0}};
CONSTANT static UpTemplate const ref3[] = {{1,2,1},{3,1,0}};
CONSTANT static UpTemplate const ref4[] = {{2,2,1},{1,1,0}};
CONSTANT static UpTemplate const ref5[] = {{3,2,1},{2,1,0}};
CONSTANT static UpTemplate const* const ref_[] = {
  ref0,ref1,ref2,ref3,ref4,ref5};
CONSTANT static UpTemplate const* const* const r_u_[] = {
  rve_,ref_};
CONSTANT static UpTemplate const* const* const* const up_templates[] = {
  0,0,f_u_,r_u_};

extern I8 const degrees[DIMS][DIMS];
extern char const* const singular_names[DIMS];
extern char const* const plural_names[DIMS];
