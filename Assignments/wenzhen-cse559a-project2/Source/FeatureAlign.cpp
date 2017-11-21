///////////////////////////////////////////////////////////////////////////
//
// NAME
//  FeatureAlign.h -- image registration using feature matching
//
// SEE ALSO
//  FeatureAlign.h      longer description
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE576 Spring 2005)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "FeatureAlign.h"
#include <math.h>

/******************* TO DO *********************
* alignPair:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*               *NOTE* Each match in 'matches' contains two feature ids of matching features, id1 (in f1) and id2 (in f2).
*               These ids are 1-based indices into the feature arrays,
*               so you access the appropriate features as f1[id1-1] and f2[id2-1].
*		m: motion model
*		f: focal length
*		nRANSAC: number of RANSAC iterations
*		RANSACthresh: RANSAC distance threshold
*		M: transformation matrix (output)
*	OUTPUT:
*		repeat for nRANSAC iterations:
*			choose a minimal set of feature matches
*			estimate the transformation implied by these matches
*			count the number of inliers
*		for the transformation with the maximum number of inliers,
*		compute the least squares motion estimate using the inliers,
*		and store it in M
*/
int alignPair(const FeatureSet &f1, const FeatureSet &f2,
              const vector<FeatureMatch> &matches, MotionModel m, float f,
              int nRANSAC, double RANSACthresh, CTransform3x3& M)
{
    // BEGIN TODO
    // write this entire method

    // inter-image transform matrix T
    CTransform3x3 T;
    T[0][0] = 1;    T[0][1] = 0;    T[0][2] = 0;
    T[1][0] = 0;    T[1][1] = 1;    T[1][2] = 0;
    T[2][0] = 0;    T[2][1] = 0;    T[2][2] = 1;

    int matchID = 0, id1 = 0, id2 = 0;
    vector<int> inliers, inliersT;
    int most = 0, idMost = 0;
    
    // repeate nRANSAC trials
    for(int i = 0; i < nRANSAC; i++){
        matchID = rand() % (matches.size());
        id1 = matches[matchID].id1;
        id2 = matches[matchID].id2;

        T[0][2] = f2[id2 - 1].x - f1[id1 - 1].x;
        T[1][2] = f2[id2 - 1].y - f1[id1 - 1].y;

        // count how many of the feature matches agree
        int num = countInliers(f1, f2, matches, m, f, T, RANSACthresh, inliersT);

        // the motion estimate with the largest number of inliers
        if(num > most){
            idMost = matchID;
            inliers = inliersT;
            most = num;
        }
    }

    leastSquaresFit(f1, f2, matches, m, f, inliers, M);
    // END TODO
    return 0;
}

/******************* TO DO *********************
* countInliers:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*               *NOTE* Each match contains two feature ids of matching features, id1 (in f1) and id2 (in f2).
*               These ids are 1-based indices into the feature arrays,
*               so you access the appropriate features as f1[id1-1] and f2[id2-1].
*		m: motion model
*		f: focal length
*		M: transformation matrix
*		RANSACthresh: RANSAC distance threshold
*		inliers: inlier feature IDs
*	OUTPUT:
*		transform the matched features in f1 by M
*
*		count the number of matching features for which the transformed
*		feature f1[id1-1] is within SSD distance RANSACthresh of its match
*		f2[id2-1]
*
*		store the indices of these matches in inliers
*
*		
*/
int countInliers(const FeatureSet &f1, const FeatureSet &f2,
                 const vector<FeatureMatch> &matches, MotionModel m, float f,
                 CTransform3x3 M, double RANSACthresh, vector<int> &inliers)
{
    inliers.clear();
    int count = 0;

    for (unsigned int i=0; i<(int) matches.size(); i++) {
        // BEGIN TODO
        // determine if the ith matched feature f1[id1-1], when transformed by M,
        // is within RANSACthresh of its match in f2
        //
        // if so, increment count and append i to inliers
        //
        // *NOTE* Each match contains two feature ids of matching features, id1 and id2.
        //        These ids are 1-based indices into the feature arrays,
        //        so you access the appropriate features as f1[id1-1] and f2[id2-1].
        int id_1 = matches[i].id1;
        int id_2 = matches[i].id2;

        Feature feature1 = f1[id_1 - 1];
        Feature feature2 = f2[id_2 - 1];

        double tX = feature1.x + M[0][2];
        double tY = feature1.y + M[1][2];

        double dx = feature2.x - tX;
        double dy = feature2.y - tY;

        double distance = sqrt(dx * dx + dy * dy);

        if(distance < RANSACthresh){
            inliers.push_back(i);
            count++;       
        }
        // END TODO
    }

    return count;
}

/******************* TO DO *********************
* leastSquaresFit:
*	INPUT:
*		f1, f2: source feature sets
*		matches: correspondences between f1 and f2
*		m: motion model
*		f: focal length
*		inliers: inlier match indices (indexes into 'matches' array)
*		M: transformation matrix (output)
*	OUTPUT:
*		compute the transformation from f1 to f2 using only the inliers
*		and return it in M
*/
int leastSquaresFit(const FeatureSet &f1, const FeatureSet &f2,
                    const vector<FeatureMatch> &matches, MotionModel m, float f,
                    const vector<int> &inliers, CTransform3x3& M)
{
    // for project 2, the transformation is a translation and
    // only has two degrees of freedom
    //
    // therefore, we simply compute the average translation vector
    // between the feature in f1 and its match in f2 for all inliers
    double u = 0;
    double v = 0;

    for (int i=0; i<inliers.size(); i++) {
        double xTrans, yTrans;

        // BEGIN TODO
        // compute the translation implied by the ith inlier match
        // and store it in (xTrans,yTrans)
        xTrans = f2[matches[inliers[i]].id2 - 1].x - f1[matches[inliers[i]].id1 - 1].x;
		yTrans = f2[matches[inliers[i]].id2 - 1].y - f1[matches[inliers[i]].id1 - 1].y;

        // END TODO
        u += xTrans;
        v += yTrans;
    }

    u /= inliers.size();
    v /= inliers.size();

    M[0][0] = 1;
    M[0][1] = 0;
    M[0][2] = -1.*u;
    M[1][0] = 0;
    M[1][1] = 1;
    M[1][2] = -1.*v;
    M[2][0] = 0;
    M[2][1] = 0;
    M[2][2] = 1;

    return 0;
}
