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
CONSTANT static I8 const rfv3[] = {0,3,2};
CONSTANT static I8 const* const rfv_[] = {
  rfv0,rfv1,rfv2,rfv3};
CONSTANT static I8 const* const* const r_v_[] = {
  0,rev_,rfv_};
CONSTANT static I8 const* const* const* const simplices[] = {
  0,0,f_v_,r_v_};
