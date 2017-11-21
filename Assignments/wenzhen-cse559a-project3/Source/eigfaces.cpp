
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: eigfaces.cpp                                                                         //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <math.h>
#include <algorithm>

using namespace std;

EigFaces::EigFaces()
:
Faces()
{
	//empty
}

EigFaces::EigFaces(int count, int width, int height)
:
Faces(count, width, height)
{
	//empty
}

void EigFaces::projectFace(const Face& face, Vector& coefficients) const
{
	if (face.getWidth()!=width || face.getHeight()!=height) {
		throw Error("Project: Face to project has different dimensions");
	}

	coefficients.resize(getSize());
	// ----------- TODO #2: compute the coefficients for the face and store in coefficients.
	int eigenFaceDimension = getSize();
	double a = 0.0;

	for(int i = 0; i < eigenFaceDimension; i++){
		a = 0.0;
		for(int y = 0; y < height; y++){
			for(int x = 0; x < width; x++){
				// a = (x - m) * u  (u is an eigenvector)
				double faceX, faceMean, eigenVecU;
				faceX = (face.pixel(x, y, 0));
				faceMean = average_face.pixel(x, y, 0);
				eigenVecU = (*this)[i].pixel(x, y, 0);
				a += (faceX - faceMean) * eigenVecU;
			}
		}
		coefficients[i] = a;
	}
}

void EigFaces::constructFace(const Vector& coefficients, Face& result) const
{	
	// ----------- TODO #3: construct a face given the coefficients
	int eigenFaceDimension = getSize();
	// x = m + sum[a_i * u_i, {i, 0, m-1}]
	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++){
			result.pixel(x, y, 0) = average_face.pixel(x, y, 0);
			for(int i = 0; i < eigenFaceDimension; i++){
				result.pixel(x, y, 0) += coefficients[i] * (*this)[i].pixel(x, y, 0);
			}
		}
	}
}

bool EigFaces::isFace(const Face& face, double max_reconstructed_mse, double& mse) const
{
	// ----------- TODO #4: Determine if an image is a face and return true if it is. Return the actual
	// MSE you calculated for the determination in mse
	// Be sure to test this method out with some face images and some non face images
	// to verify it is working correctly.

	Face X(width, height);
	int dimension = width * height;
	double squareErrorPerPixel;
	mse = 0.0;
	Vector coefficients;

	projectFace(face, coefficients);
	constructFace(coefficients, X);

	for(int y = 0; y < height; y++){
		for(int x = 0; x < width; x++){
			squareErrorPerPixel = face.pixel(x, y, 0) - X.pixel(x, y, 0);
			mse += squareErrorPerPixel * squareErrorPerPixel;
		}
	}
	mse = mse / dimension;

	if(mse < max_reconstructed_mse){
		return true;
	} else{
		return false;
	}
}

bool EigFaces::verifyFace(const Face& face, const Vector& user_coefficients, double max_coefficients_mse, double& mse) const
{
	// ----------- TODO #5 : Determine if face is the same user give the user's coefficients.
	// return the MSE you calculated for the determination in mse.
	Vector coefficients;
	mse = 0.;
	int eigenFaceDimension = getSize();

	projectFace(face, coefficients);
	for (int i = 0; i < eigenFaceDimension; i++){
		mse += pow((user_coefficients[i] - coefficients[i]), 2);
	}
	mse /= eigenFaceDimension;
	if(mse <= max_coefficients_mse){
		return true;
	}else{
		return false;
	}
}

void EigFaces::recognizeFace(const Face& face, Users& users) const
{
	// ----------- TODO #6: Sort the users by closeness of match to the face
	Vector coefficients;

	double mse = 0.;
	int eigenFaceDimension = getSize();

	projectFace(face, coefficients);
	for(int i = 0; i < users.getSize(); i++){
		mse = 0.;
		for(int j = 0; j < eigenFaceDimension; j++){
			mse += pow((users[i][j] - coefficients[j]), 2);
		}
		mse /= eigenFaceDimension;
		users[i].setMse(mse);
	}
	users.sort();
}

bool sortFacePosition(FacePosition fp1, FacePosition fp2){
	return fp1.error < fp2.error;
}

void EigFaces::findFace(const Image& img, double min_scale, double max_scale, double step, int n, bool crop, Image& result) const
{
	// ----------- TODO #7: Find the faces in Image. Search image scales from min_scale to max_scale inclusive,
	// stepping by step in between. Find the best n faces that do not overlap each other. If crop is true,
	// n is one and you should return the cropped original img in result. The result must be identical
	// to the original besides being cropped. It cannot be scaled and it must be full color. If crop is
	// false, draw green boxes (use r=100, g=255, b=100) around the n faces found. The result must be
	// identical to the original image except for the addition of the boxes.

	vector<FacePosition> bestMatches (n); // best n faces
	FacePosition bestMatch;
	vector<FacePosition>::iterator it;

	int cropWidth = width;
	int cropHeight = height;

	int w = img.getWidth();
	int h = img.getHeight();

	Image scaledImage;
	Face subface (width, height);

	double max_mse = DBL_MAX;
	double threshold = sqrt(width * width + height * height);
	// scale the face in the outermost loop
	for(double scale = min_scale; scale < max_scale; scale += step){
		cout << "scale = " << scale << endl;
		// scale the face
		scaledImage.resize((int) w * scale, (int) h * scale);
		// sample
		img.resample(scaledImage);

		int scaledImageWidth = scaledImage.getWidth();
		int scaledImageHeight = scaledImage.getHeight();

		for(int y = 0; y < scaledImageHeight; y++){
			for(int x = 0; x < scaledImageWidth; x++){

				subface.subimage(x, x + subface.getWidth() - 1, y, y + subface.getHeight() - 1, scaledImage, false);
				double mse;
				bool isface = isFace(subface, max_mse, mse);

				if(isface){
					FacePosition tempPos = FacePosition();
					tempPos.x = x;
					tempPos.y = y;
					tempPos.scale = scale;
					tempPos.error = mse;

					it = bestMatches.begin();
					if(it == bestMatches.end()){
						bestMatches.push_back(tempPos);
					}
					else{
						bool overlap = false;
						for(vector<FacePosition>::iterator it = bestMatches.begin(); it != bestMatches.end(); it++){
							double xdistance = ((*it).x/(*it).scale) - (tempPos.x/tempPos.scale);
							double ydistance = ((*it).y/(*it).scale) - (tempPos.y/tempPos.scale);
							double distance = sqrt(xdistance * xdistance +ydistance * ydistance);
							// overlap
							if (distance < threshold){
								overlap = true;
								if((*it).error > tempPos.error){
									*it = tempPos;
									sort(bestMatches.begin(), bestMatches.end(), sortFacePosition);
								}
								break;
							}
						}
						// overlap = false:
						if(!overlap){
							bestMatches.push_back(tempPos);
							sort(bestMatches.begin(), bestMatches.end(), sortFacePosition);
							// bestMatches.pop_back(); // delete the last one
						}
					}
				}//end of the dealing with subface like a face
			}
		}//end of loop for scaled image
	}// specific scale

	while(bestMatches.size() > n){
		bestMatches.pop_back();
	}

	for(vector<FacePosition>::iterator it = bestMatches.begin(); it != bestMatches.end(); it++){
		cout << (*it).error << " " ;
	}
	cout << endl;



	if(!crop){
		result.resize(w, h, img.getColors());
		img.resample(result);
		int x1, y1, x2, y2;
		for(vector<FacePosition>:: iterator it = bestMatches.begin(); it != bestMatches.end(); it++){
			x1 = (int)((*it).x/(*it).scale);
			y1 = (int)((*it).y/(*it).scale);
			x2 = (int)(((*it).x+cropWidth)/((*it).scale));
			y2 = (int)(((*it).y+cropHeight)/((*it).scale));

			result.line(x1, y1, x2, y1, 100, 256, 100);
			result.line(x2, y1, x2, y2, 100, 256, 100);
			result.line(x1, y1, x1, y2, 100, 256, 100);
			result.line(x1, y2, x2, y2, 100, 256, 100);
		}
	}
	else{
		FacePosition win = bestMatches.front();
		img.crop(win.x/win.scale, win.y/win.scale, (win.x + cropWidth)/win.scale, (win.y+cropHeight)/win.scale, result);
	}
}



void EigFaces::morphFaces(const Face& face1, const Face& face2, double distance, Face& result) const
{
	// TODO (extra credit): MORPH along *distance* fraction of the vector from face1 to face2 by
	// interpolating between the coefficients for the two faces and reconstructing the result.
	// For example, distance 0.0 will approximate the first, while distance 1.0 will approximate the second.
	// Negative distances are ok two.

}

const Face& EigFaces::getAverage() const
{
	return average_face;
}

void EigFaces::setAverage(const Face& average)
{
	average_face=average;
}



