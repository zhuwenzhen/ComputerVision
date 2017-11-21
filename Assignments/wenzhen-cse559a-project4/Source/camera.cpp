/***************************************************************
 * CS4670/5670, Fall 2012 Project 4
 * File to be modified #4:
 * camera.cpp
 *	     - routines for finding the camera position and projection matrix
 *       - routine for inverting the scene
 **************************************************************/
#pragma warning(disable : 4996)

#include "ImgView.h"

#include "svmCMath.h"
using namespace std;

//
// TODO 9: computeCameraParameters()
// Compute the camera position
void ImgView::computeCameraParameters()
{
  if (refPointOffPlane == NULL) {
    fl_message("Reference height must be set\n");
    return;
  }

  if (!homographyComputed) {
    fl_message("Must first compute reference homography\n");
    return;
  }

  // Compute the position of the camera.  You'll do this in two steps,
  // first finding the height of the camera, then finding the x and y coordinates
  // of the camera, as described in the project page.

  /******** BEGIN TODO Part 1 ********/
  // Compute the height of the camera, store in z_cam, and the x and y coordinates of the camera,
  // storing in x_cam and y_cam
  double z_cam = 0.0;
  double x_cam = 0.0, y_cam = 0.0;

  // height of camera 
  Vec3d vx=Vec3d(xVanish.u,xVanish.v,xVanish.w);
  Vec3d vy=Vec3d(yVanish.u,yVanish.v,yVanish.w);
  Vec3d horizonV = cross(vx,vy);
  // vertical line
  Mat3d matH=Mat3d(
        H[0][0],H[0][1],H[0][2],
        H[1][0],H[1][1],H[1][2],
        H[2][0],H[2][1],H[2][2]);
  Vec3d r = Vec3d(refPointOffPlane->u,refPointOffPlane->v,refPointOffPlane->w);
  Vec3d pos3dOfRef = Vec3d(refPointOffPlane->X,refPointOffPlane->Y,1);
  Vec3d r0 = matH * pos3dOfRef;

  // find the intersection
  Vec3d vertical = cross(r,r0);
  Vec3d height = cross(horizonV, vertical);

  SVMPoint hp,rp;

  // normalized
  r0 /= r0[2];
  height /= height[2];

  hp.u = height[0];
  hp.v = height[1];
  hp.w = height[2];

  rp.u = r0[0];
  rp.v = r0[1];
  rp.w = r0[2];
  rp.X = refPointOffPlane->X;
  rp.Y = refPointOffPlane->Y;
  rp.Z = 0;
  rp.W = 1;

  rp.known(true);

  pntSelStack.push_back(&rp);
  pntSelStack.push_back(&hp);
  sameXY();
  pntSelStack.pop_back();
  pntSelStack.pop_back();

  z_cam = hp.Z;

  // camera point
  Vec3d vz=Vec3d(zVanish.u,zVanish.v,zVanish.w);
  vz /= vz[2];

  SVMPoint cameraPt;
  cameraPt.u = vz[0];
  cameraPt.v = vz[1];
  cameraPt.w = 1;

  pntSelStack.push_back(&rp);
  pntSelStack.push_back(&cameraPt);
  sameZPlane();
  pntSelStack.pop_back();
  pntSelStack.pop_back();

  x_cam = cameraPt.X;
  y_cam = cameraPt.Y;

 /******** END TODO Part 1 ********/
 camPos[0] = x_cam;
 camPos[1] = y_cam;
 camPos[2] = z_cam;

 printf("Camera is at [ %0.3f %0.3f %0.3f ]\n", camPos[0], camPos[1], camPos[2]);

 camComputed = true;
}

//
// EXTRA CREDIT 10: invertScene(double zScale)
//   "invert" the scene, creating the reverse perspective
//   You'll need to compute the right 4x4 transformation matrix for inverting the
//     scene (see webpage for more details)
//   zScale is the factor by which you'll scale the z-direction of the inverted scene
//     (see webpage for more details)
void ImgView::invertScene(double zScale)
{

  fl_message("Extra credit: invert scene\n");
  printf("Extra credit: invert scene\n");

  if (sceneInverted) {
    fl_message("Scene has already been inverted\n");
    return;
  }

  if (!camComputed) {
    fl_message("Compute projection matrix first\n");
    return;
  }

  /******** BEGIN TODO ********/
  // Compute the transform which will "invert" the scene, creating a reverse perspective
  // Store the results in the 4x4 matrix 'inv'.  You'll need the camera translation and rotation,
  // which are stored in campos and camR, respectively
  Mat4d inv;

  Mat4d L = Mat4d(1, 0, 0, 0, 
                  0, 1, 0, 0, 
                  0, 0, 0, zScale, 
                  0, 0, 1, 0);
  Mat4d M = Mat4d(camR[0][0],camR[0][1],camR[0][2],0,
	                camR[1][0],camR[1][1],camR[1][2],0,
				          camR[2][0],camR[2][1],camR[2][2],0,
				          0,0,0,1);

  // cout << camR[0][0]  << " ";
  // cout << camR[0][1]  << " ";
  // cout << camR[0][2]  << endl;
  // cout << camR[1][0]  << " ";
  // cout << camR[1][1]  << " ";
  // cout << camR[1][2]  << endl;
  // cout << camR[2][0]  << " ";
  // cout << camR[2][1]  << " ";
  // cout << camR[2][2]  << endl;

  Mat4d R = Mat4d(1,0,0,-1*camPos[0],
	                0,1,0,-1*camPos[1],
				          0,0,1,-1*camPos[2],
				          0,0,0,1);

  // cout << camPos[0]  << endl;
  // cout << camPos[1]  << endl;
  // cout << camPos[2]  << endl;

  inv = L * M * R;

 /********* END TODO ********/


 // This loop applies the matrix to each point in the scene
 int index = 0;
 CTypedPtrDblElement <SVMPoint> *pntNode = pntList.GetHeadPtr();
 while (!pntList.IsSentinel(pntNode))
   {
     SVMPoint *pnt = pntNode->Data();
     index++;

     if (!pnt->known()) {
       pntNode = pntNode->Next();
       continue;
     }

     Vec4d p(pnt->X, pnt->Y, pnt->Z, pnt->W);
     Vec4d pNew = inv * p;

     printf("inverting point %d:  %0.3f %0.3f %0.3f => %0.3f %0.3f %0.3f\n",
            index, pnt->X, pnt->Y, pnt->Z, pNew[0] / pNew[3], pNew[1] / pNew[3], pNew[2] / pNew[3]);

     pnt->X = pNew[0] / pNew[3];
     pnt->Y = pNew[1] / pNew[3];
     pnt->Z = pNew[2] / pNew[3];
     pnt->W = 1.0;

     pntNode = pntNode->Next();
   }

 CTypedPtrDblElement <SVMPolygon> *plyNode = plyList.GetHeadPtr();
 // Force populate homography
 while (!plyList.IsSentinel(plyNode)) {
   SVMPolygon *ply = plyNode->Data();
   populateHomography(*ply, NULL);
   ply->isHomographyPopulated = true;
   plyNode=plyNode->Next();
 }

 sceneInverted = true;
}






/* Use a 180 rotation to fix up the intrinsic matrix */
void FixIntrinsics(double *P, double *K, double *R, double *t)
{
#if 0
  /* Check the parity along the diagonal */
  int neg = (K[0] < 0.0) + (K[4] < 0.0) + (K[8] < 0.0);

  /* If odd parity, negate the instrinsic matrix */
  if ((neg % 2) == 1) {
    matrix_scale(3, 3, K, -1.0, K);
    matrix_scale(3, 4, P, -1.0, P);
  }

  /* Now deal with case of even parity */
  double fix[9];
  matrix_ident(3, fix);
  double tmp[9], tmp2[12];

  if (K[0] < 0.0 && K[4] < 0.0) {
    fix[0] = -1.0;
    fix[4] = -1.0;
  } else if (K[0] < 0.0) {
    fix[0] = -1.0;
    fix[8] = -1.0;
  } else if (K[4] < 0.0) {
    fix[4] = -1.0;
    fix[8] = -1.0;
  } else {
    /* No change needed */
  }

  matrix_product(3, 3, 3, 3, K, fix, tmp);
  memcpy(K, tmp, sizeof(double) * 3 * 3);

  double Kinv[9];
  matrix_invert(3, K, Kinv);

  matrix_product(3, 3, 3, 4, Kinv, P, tmp2);

  memcpy(R + 0, tmp2 + 0, sizeof(double) * 3);
  memcpy(R + 3, tmp2 + 4, sizeof(double) * 3);
  memcpy(R + 6, tmp2 + 8, sizeof(double) * 3);

  t[0] = tmp2[3];
  t[1] = tmp2[7];
  t[2] = tmp2[11];
#endif
}


void ImgView::decomposePMatrix() {
#if 0
  double KRinit[9], Kinit[9], Rinit[9], tinit[3];
  memcpy(KRinit + 0, camP[0], 3 * sizeof(double));
  memcpy(KRinit + 3, camP[1], 3 * sizeof(double));
  memcpy(KRinit + 6, camP[2], 3 * sizeof(double));

  dgerqf_driver(3, 3, KRinit, Kinit, Rinit);

  /* We want our intrinsics to have a certain form */
  FixIntrinsics(camP[0], Kinit, Rinit, tinit);
  matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

  memcpy(camK[0], Kinit, 9 * sizeof(double));

  camR = Mat4d(1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0,
               0.0, 0.0, 1.0, 0.0,
               0.0, 0.0, 0.0, 1.0);

  memcpy(camR[0], Rinit + 0, 3 * sizeof(double));
  memcpy(camR[1], Rinit + 3, 3 * sizeof(double));
  memcpy(camR[2], Rinit + 6, 3 * sizeof(double));

  printf("R:\n"
         "[ %0.3f %0.3f %0.3f\n"
         "  %0.3f %0.3f %0.3f\n"
         "  %0.3f %0.3f %0.3f ]\n",
         camR[0][0], camR[0][1], camR[0][2],
         camR[1][0], camR[1][1], camR[1][2],
         camR[2][0], camR[2][1], camR[2][2]);

  printf("K:\n"
         "[ %0.3f %0.3f %0.3f\n"
         "  %0.3f %0.3f %0.3f\n"
         "  %0.3f %0.3f %0.3f ]\n",
         camK[0][0], camK[0][1], camK[0][2],
         camK[1][0], camK[1][1], camK[1][2],
         camK[2][0], camK[2][1], camK[2][2]);
#endif
}
