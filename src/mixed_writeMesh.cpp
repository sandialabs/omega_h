#include "MeshSim.h"
#include "SimDiscrete.h"
#include "SimMessages.h"
#include "SimError.h"
#include "SimErrorCodes.h"
#include "SimMeshingErrorCodes.h"
#include "SimDiscreteErrorCodes.h"
#include <iostream>
using namespace std;

#include "Omega_h_comm.hpp"
#include "Omega_h_mesh.hpp"
#include "Omega_h_file.hpp"
using namespace Omega_h;

void messageHandler(int type, const char *msg)
{
  switch (type) {
  case Sim_InfoMsg:
    cout<<"Info: "<<msg<<endl;
    break;
  case Sim_DebugMsg:
    cout<<"Debug: "<<msg<<endl;
    break;
  case Sim_WarningMsg:
    cout<<"Warning: "<<msg<<endl;
    break;
  case Sim_ErrorMsg:
    cout<<"Error: "<<msg<<endl;
    break;
  }
  return;
}

void get_pyramid(int *numVerts, int *numElems,
                 int *elementType, int *elementData, double *coords)
{
    numVerts = new int[1];
    *numVerts = 5;
    numElems = new int[1];
    *numElems = 1;
    elementType = new int[numElems[0]];
    *elementType = 11;
    elementData = new int[numElems[0]*5];
    elementData[0] = 0; elementData[1] = 1; elementData[2] = 2;
    elementData[3] = 3; elementData[4] = 4;
    coords = new double[numVerts[0]*3];
    coords = new double[numVerts[0]*3];
    coords[0] = 0; coords[1] = 0; coords[2] = 1;
    coords[3] = 1; coords[4] = 0; coords[5] = 1;
    coords[6] = 1; coords[7] = 1; coords[8] = 1;
    coords[9] = 0; coords[10]= 1; coords[11]= 1;
    coords[12]=0.5;coords[13]=0.5;coords[14]= 2;
}

int main(int argc, char *argv[]) {
  auto lib = Library(&argc, &argv);
  auto comm = lib.world();
  auto mesh = Mesh(comm->library());
  try {
    Sim_logOn("importData1.log");
    MS_init();
    Sim_readLicenseFile(0);

    pMesh meshtest;

/*
    // For tet on wedge, plus hex and pyramid adjacenct to wedge
    int numVerts = 12;
    double coords[12*3] = {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000, 
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                           1.000, 1.500, 0.500,
                           0.000, -1.000, 0.000,
                           1.000, -1.000, 0.000,
                           0.000, -1.000, 1.000,
                           1.000, -1.000, 1.000
                          };
    int numElems = 4;
    int elementType[4] = {12, 10, 11, 13};
    int elementData[6+4+5+8] = {0,1,2,3,4,5, 3,4,5,6, 5,4,1,2,7,
8,9,1,0,10,11,4,3};
    pVertex vReturn[12]; // array of created vertices
    pEntity eReturn[4]; // array of created entities
    //
*/
/*
    //For pyramid on hex
    int numVerts = 9;
    double coords[9*3] =  {
                           0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           1.000, 1.000, 0.000,
                           0.000, 1.000, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000,
                           0.500, 0.500, 2.000
                          };
    int numElems = 2;
    int elementType[2] = {13, 11};
    int elementData[8+5] = {0,1,2,3,4,5,6,7, 4,5,6,7,8};
    pVertex vReturn[9];
    pEntity eReturn[2];
    //For pyramid
    int numVerts = 5;
    double coords[5*3] =  {
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000,
                           0.500, 0.500, 2.000
                          };
    int numElems = 1;
    int elementType[1] = {11};
    int elementData[5] = {0,1,2,3,4};
    pVertex vReturn[5];
    pEntity eReturn[1];
*/

    int numVerts, numElems, elementType, elementData;
    double coords;
    get_pyramid(&numVerts, &numElems, &elementType, &elementData, &coords);
    pVertex vReturn[5];
    pEntity eReturn[1];
/*
    //For single hex
    int numVerts = 8;
    double coords[8*3] =  {
                           0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           1.000, 1.000, 0.000,
                           0.000, 1.000, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           1.000, 1.000, 1.000,
                           0.000, 1.000, 1.000
                          };
    int numElems = 1;
    int elementType[1] = {13};
    int elementData[8+5] = {0,1,2,3,4,5,6,7};
    pVertex vReturn[8];
    pEntity eReturn[1];
*/
/*
    //For tet on wedge
    int numVerts = 7;
    double coords[7*3] =  {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                           0.500, 0.289, 1.866,
                          };
    int numElems = 2;
    int elementType[2] = {12, 10};
    int elementData[6+4] = {0,1,2,3,4,5,3,4,5,6};
    pVertex vReturn[7];
    pEntity eReturn[2];
*/
/*
    //For single wedge
    int numVerts = 6;
    double coords[6*3] =  {0.000, 0.000, 0.000,
                           1.000, 0.000, 0.000,
                           0.500, 0.866, 0.000,
                           0.000, 0.000, 1.000,
                           1.000, 0.000, 1.000,
                           0.500, 0.866, 1.000,
                          };
    int numElems = 1;
    int elementType[1] = {12};
    int elementData[6] = {0,1,2,3,4,5};
    pVertex vReturn[6];
    pEntity eReturn[1];
*/
    SimDiscrete_start(0);
    Sim_setMessageHandler(messageHandler);
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    meshtest = M_new(0,0);
    if(M_importFromData(meshtest,numVerts,&coords,numElems,
      &elementType,&elementData,vReturn,eReturn,progress)) {
      cerr<<"Error importing mesh data"<<endl;
      M_release(meshtest);
    }

    pDiscreteModel modeltest = DM_createFromMesh(meshtest, 0, progress);
    if(!modeltest) {
      cerr<<"Error creating Discrete model from mesh"<<endl;
      M_release(meshtest);
    }

    DM_findEdgesByFaceNormals(modeltest, 0, progress);
    DM_eliminateDanglingEdges(modeltest, progress);
    if(DM_completeTopology(modeltest, progress)) {
      cerr<<"Error completing Discrete model topology"<<endl;
      M_release(meshtest);
      GM_release(modeltest);
    }

    GM_write(modeltest,"/users/joshia5/simmodeler/Example_pym.smd",0,progress);
    M_write(meshtest,"/users/joshia5/simmodeler/Example_pym.sms", 0,progress);

    M_release(meshtest);
    GM_release(modeltest);
    Progress_delete(progress);
    SimDiscrete_stop(0);
    Sim_unregisterAllKeys();
    MS_exit();
    Sim_logOff();

  } catch (pSimError err) {
    cerr<<"SimModSuite error caught:"<<endl;
    cerr<<"  Error code: "<<SimError_code(err)<<endl;
    cerr<<"  Error string: "<<SimError_toString(err)<<endl;
    SimError_delete(err);
  } catch (...) {
    cerr<<"Unhandled exception caught"<<endl;
  }
  return 0;
}
