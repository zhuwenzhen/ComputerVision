#include <assert.h>
#include <math.h>
#include <FL/Fl.H>
#include <FL/Fl_Image.H>
#include "features.h"
#include "ImageLib/FileIO.h"               
#include <algorithm>
#include <iterator>
#include <numeric>
#include <functional>

#define PI 3.14159265358979323846

// Compute features of an image.
bool computeFeatures(CFloatImage &image, FeatureSet &features, int featureType) {
	// TODO: Instead of calling dummyComputeFeatures, write your own
	// feature computation routines and call them here.
	switch (featureType) {
	case 1:
		cout << "Simple Feature" << endl;
		ComputeSimpleFeatures(image, features);
		break;
	case 2:
		cout << "Harris Features" << endl;
		ComputeHarrisFeatures(image, features);		
		break;
	case 3:
		cout << "my Feature" << endl;
		ComputeMyFeatures(image, features);
		break;
		
	default:
		return false;
	}

	// This is just to make sure the IDs are assigned in order, because
	// the ID gets used to index into the feature array.
	for (unsigned int i = 0; i<features.size(); i++) {
		features[i].id = i + 1;
	}
	return true;
}

// Perform a query on the database.  This simply runs matchFeatures on
// each image in the database, and returns the feature set of the best
// matching image.
bool performQuery(const FeatureSet &f, const ImageDatabase &db, int &bestIndex, vector<FeatureMatch> &bestMatches, double &bestScore, int matchType) {
	// Here's a nice low number.
	bestScore = -1e100;

	vector<FeatureMatch> tempMatches;
	double tempScore;

	for (unsigned int i = 0; i<db.size(); i++) {
		if (!matchFeatures(f, db[i].features, tempMatches, tempScore, matchType)) {
			return false;
		}

		if (tempScore > bestScore) {
			bestIndex = i;
			bestScore = tempScore;
			bestMatches = tempMatches;
		}
	}

	return true;
}

// Match one feature set with another.
bool matchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore, int matchType) {
	// TODO: We have given you the ssd matching function, you must write your own
	// feature matching function for the ratio test.

	printf("\nMatching features.......\n");

	switch (matchType) {
	case 1:
		ssdMatchFeatures(f1, f2, matches, totalScore);
		return true;
	case 2:
		ratioMatchFeatures(f1, f2, matches, totalScore);
		return true;
	default:
		return false;
	}
}

// Evaluate a match using a ground truth homography.  This computes the
// average SSD distance between the matched feature points and
// the actual transformed positions.
double evaluateMatch(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9]) {
	double d = 0;
	int n = 0;

	double xNew;
	double yNew;

	unsigned int num_matches = matches.size();
	for (unsigned int i = 0; i<num_matches; i++) {
		int id1 = matches[i].id1;
		int id2 = matches[i].id2;
		applyHomography(f1[id1 - 1].x, f1[id1 - 1].y, xNew, yNew, h);
		d += sqrt(pow(xNew - f2[id2 - 1].x, 2) + pow(yNew - f2[id2 - 1].y, 2));
		n++;
	}

	return d / n;
}

void addRocData(const FeatureSet &f1, const FeatureSet &f2, const vector<FeatureMatch> &matches, double h[9], vector<bool> &isMatch, double threshold, double &maxD) {
	double d = 0;

	double xNew;
	double yNew;

	unsigned int num_matches = matches.size();
	for (unsigned int i = 0; i<num_matches; i++) {
		int id1 = matches[i].id1;
		int id2 = matches[i].id2;
		applyHomography(f1[id1 - 1].x, f1[id1 - 1].y, xNew, yNew, h);

		// Ignore unmatched points.  There might be a better way to
		// handle this.
		d = sqrt(pow(xNew - f2[id2 - 1].x, 2) + pow(yNew - f2[id2 - 1].y, 2));
		if (d <= threshold)
		{
			isMatch.push_back(1);
		}
		else
		{
			isMatch.push_back(0);
		}

		if (matches[i].score>maxD)
			maxD = matches[i].score;
	}
}

vector<ROCPoint> computeRocCurve(vector<FeatureMatch> &matches, vector<bool> &isMatch, vector<double> &thresholds)
{
	vector<ROCPoint> dataPoints;

	for (int i = 0; i < (int)thresholds.size(); i++)
	{
		//printf("Checking threshold: %lf.\r\n",thresholds[i]);
		int tp = 0;
		int actualCorrect = 0;
		int fp = 0;
		int actualError = 0;
		int total = 0;

		int num_matches = (int)matches.size();
		for (int j = 0; j < num_matches; j++)
		{
			if (isMatch[j])
			{
				actualCorrect++;
				if (matches[j].score<thresholds[i])
				{
					tp++;
				}
			}
			else
			{
				actualError++;
				if (matches[j].score<thresholds[i])
				{
					fp++;
				}
			}
			total++;
		}

		ROCPoint newPoint;
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		newPoint.trueRate = (double(tp) / actualCorrect);
		newPoint.falseRate = (double(fp) / actualError);
		//printf("newPoints: %lf,%lf",newPoint.trueRate,newPoint.falseRate);
		dataPoints.push_back(newPoint);
	}
	return dataPoints;
}
 

// Compute silly example features.  This doesn't do anything
// meaningful.
void dummyComputeFeatures(CFloatImage &image, FeatureSet &features) {
	CShape sh = image.Shape();
	Feature f;

	for (int y = 0; y<sh.height; y++) {
		for (int x = 0; x<sh.width; x++) {
			double r = image.Pixel(x, y, 0);
			double g = image.Pixel(x, y, 1);
			double b = image.Pixel(x, y, 2);

			if ((int)(255 * (r + g + b) + 0.5) % 100 == 1) {
				// If the pixel satisfies this meaningless criterion,
				// make it a feature.
				f.type = 1;
				f.id += 1;
				f.x = x;
				f.y = y;
	
				f.data.resize(1);
				f.data[0] = r + g + b;
				features.push_back(f);
			}
		}
	}
}

void ComputeHarrisFeatures(CFloatImage &image, FeatureSet &features)
{
	//Create grayscale image used for Harris detection
	CFloatImage grayImage = ConvertToGray(image);

	//Create image to store Harris values
	CFloatImage harrisImage(image.Shape().width, image.Shape().height, 1);

	//Create image to store local maximum harris values as 1, other pixels 0
	CByteImage harrisMaxImage(image.Shape().width, image.Shape().height, 1);

	//compute Harris values puts harris values at each pixel position in harrisImage. 
	//You'll need to implement this function.
	computeHarrisValues(grayImage, harrisImage);

	// Threshold the harris image and compute local maxima.  You'll need to implement this function.
	computeLocalMaxima(harrisImage, harrisMaxImage);

	// Prints out the harris image for debugging purposes
	CByteImage tmp(harrisImage.Shape());
	convertToByteImage(harrisImage, tmp);
	WriteFile(tmp, "harris.tga");

	CFloatImage orientation(harrisMaxImage.Shape().width, harrisMaxImage.Shape().height, 1);
	computeOrientation(grayImage, orientation);

	// TO DO--------------------------------------------------------------------
	//Loop through feature points in harrisMaxImage and create feature descriptor 
	//for each point above a threshold

	for (int y = 0; y<harrisMaxImage.Shape().height; y++) {
		for (int x = 0; x<harrisMaxImage.Shape().width; x++) {
			// Skip over non-maxima
			if (harrisMaxImage.Pixel(x, y, 0) == 0)
				continue;
			//TO DO---------------------------------------------------------------------
			// Fill in feature with descriptor data here. 
			Feature f;
			f.type = 2;
			f.id += 1;
			f.x = x;
			f.y = y;
		    f.angleRadians = orientation.Pixel(x, y, 0);
			computeSimpleDescriptors(harrisImage, f);
			features.push_back(f);
		}
	}
}

// ComputeSimpleDescriptors for 5*5 window
void computeSimpleDescriptors(CFloatImage &srcImage, Feature &f) {
	int x = f.x;
	int y = f.y;

	for (int j = y - 2; j < y + 3; j++) {	
		for (int i = x - 2; i < x + 3; i++) {
			if (srcImage.Shape().InBounds(i, j)) {
				f.data.push_back(srcImage.Pixel(i, j, 0));
			}
			else {
				f.data.push_back(0.);
			}
		}
	}
}
void ComputeMyFeatures(CFloatImage &image, FeatureSet& features){
	// We only need intensity, hence convert rgb image to gray scale 
	CFloatImage grayImage = ConvertToGray(image);
	CShape sh = grayImage.Shape();

	// loop through every pixel of the grayscale image
	// make a 5*5 window, if the convolution window out of boundry, just push 0 to our feature's data

	for (int y = 0; y < sh.height; y++) {
		for (int x = 0; x < sh.width; x++) {
			Feature f;			
			vector<float> v;

			for (int j = y - 2; j < y + 3; j++) {
				for (int i = x - 2; i < x + 3; i++) {
					if (sh.InBounds(i, j)) {	
						v.push_back(grayImage.Pixel(i, j, 0));
					}
					else {
						v.push_back(0.);
					}
				}
			}

			double sum = accumulate(v.begin(), v.end(), 0.0);
			double mean = sum / v.size();

			vector<double> diff(v.size());
			transform(v.begin(), v.end(), diff.begin(),bind2nd(minus<double>(), mean));
			double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
			double stdev = sqrt(sq_sum / v.size());

			f.data.push_back(mean);
			f.data.push_back(stdev);

			f.type = 3;
			f.id += 1;
			f.x = x;
			f.y = y;
			features.push_back(f);
		}
	}
}

void ComputeSimpleFeatures(CFloatImage &image, FeatureSet& features) {
	// We only need intensity, hence convert rgb image to gray scale 
	CFloatImage grayImage = ConvertToGray(image);
	CShape sh = grayImage.Shape();

	// loop through every pixel of the grayscale image
	// make a 5*5 window, if the convolution window out of boundry, just push 0 to our feature's data

	for (int y = 0; y < sh.height; y++) {
		for (int x = 0; x < sh.width; x++) {
			Feature f;

			for (int j = y - 2; j < y + 3; j++) {
				for (int i = x - 2; i < x + 3; i++) {
					if (sh.InBounds(i, j)) {
						f.data.push_back(grayImage.Pixel(i, j, 0));
					}
					else {
						f.data.push_back(0.);
					}
				}
			}
			
			f.type = 3;
			f.id += 1;
			f.x = x;
			f.y = y;
			features.push_back(f);
		}
	}

}

void computeOrientation(CFloatImage &srcImage, CFloatImage &orientation) {
	int w = srcImage.Shape().width;
	int h = srcImage.Shape().height;

	CFloatImage dx(srcImage.Shape());
	CFloatImage dy(srcImage.Shape());

	Convolve(srcImage, dx, ConvolveKernel_SobelX);
	Convolve(srcImage, dy, ConvolveKernel_SobelY);

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			orientation.Pixel(x, y, 0) = atan2(dy.Pixel(x, y, 0), dx.Pixel(x, y, 0));
		}
	}
}


//TO DO---------------------------------------------------------------------
//Loop through the image to compute the harris corner values as described in class
// srcImage:  grayscale of original image
// harrisImage:  populate the harris values per pixel in this image
void computeHarrisValues(CFloatImage &srcImage, CFloatImage &harrisImage)
{
	int w = srcImage.Shape().width;
	int h = srcImage.Shape().height;

	CFloatImage dx(srcImage.Shape());
	CFloatImage dy(srcImage.Shape());

	CFloatImage dxx(srcImage.Shape());
	CFloatImage dyy(srcImage.Shape());
	CFloatImage dxy(srcImage.Shape());

	Convolve(srcImage, dx, ConvolveKernel_SobelX);
	Convolve(srcImage, dy, ConvolveKernel_SobelY);

	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			dxx.Pixel(x, y, 0) = dx.Pixel(x, y, 0)*dx.Pixel(x, y, 0);
			dyy.Pixel(x, y, 0) = dy.Pixel(x, y, 0)*dy.Pixel(x, y, 0);
			dxy.Pixel(x, y, 0) = dx.Pixel(x, y, 0)*dy.Pixel(x, y, 0);
		}
	}

	CFloatImage weight(5, 5, 1);
	for (int y = 0; y < 5; y++)
	{
		for (int x = 0; x < 5; x++)
		{
			weight.Pixel(x, y, 0) = gaussian5x5[y*5 + x];
		}
	}
	
	CFloatImage H11(srcImage.Shape());
	CFloatImage H22(srcImage.Shape());
	CFloatImage H12(srcImage.Shape());

	Convolve(dxx, H11, weight);
	Convolve(dyy, H22, weight);
	Convolve(dxy, H12, weight);

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			// TODO:  Compute the harris score for 'srcImage' at this pixel and store in 'harrisImage'.  See the project
			//   page for pointers on how to do this
			double h11 = H11.Pixel(x, y, 0);
			double h22 = H22.Pixel(x, y, 0);
			double h12 = H12.Pixel(x, y, 0);
			if (h11 + h22 == 0) {
				harrisImage.Pixel(x, y, 0) = 0.;
			}
			else {
				float f = (h11 * h22 - h12 * h12) / (h11 + h22);
				harrisImage.Pixel(x, y, 0) = f;
			}
		}
	}
}

// TO DO---------------------------------------------------------------------
// Loop through the harrisImage to threshold and compute the local maxima in a neighborhood
// srcImage:  image with Harris values
// destImage: Assign 1 to a pixel if it is above a threshold and is the local maximum in 3x3 window, 0 otherwise.
//    You'll need to find a good threshold to use.
void computeLocalMaxima(CFloatImage &srcImage, CByteImage &destImage)
{
	float threshold = 0.5;

	int w = srcImage.Shape().width;
	int h = srcImage.Shape().height;

	vector<float> window;

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			// compute the local maxima, use window to store 5*5 = 25 values of pixel intensity
			vector<float> window;
			for (int j = y - 2; j < y + 3; j++) {
				for (int i = x - 2; i < x + 3; i++) {
					float value;
					if (srcImage.Shape().InBounds(i, j)) {
						value = srcImage.Pixel(i, j, 0);
					}
					else {
						value = 0;
					}
					window.push_back(value);
				}
			}
			// use max_element to get the maximum 
			double local_max = *max_element(window.begin(), window.end());
			/*for (vector<float>::iterator iter = window.begin(); iter != window.end(); iter++) {
				cout << *iter << " ";
			}
			cout << endl;
			cout << "local max " << local_max << endl;*/
			// write values to the destImage
			for (int j = y - 2; j < y + 3; j++) {
				for (int i = x - 2; i < x + 3; i++) {				
					// Assign 1 to a pixel if it is above a threshold and is the local maximum
					if (srcImage.Shape().InBounds(i, j)) {
						float p_val = srcImage.Pixel(i, j, 0);					
						//	cout << "(" << x << ", " << y << "), " << p_val << " " << local_max << endl;
						if (p_val == local_max && p_val > threshold) {
							destImage.Pixel(i, j, 0) = 1;
						}
						else {
							destImage.Pixel(i, j, 0) = 0;
						}
					}
					else {
						continue;
					}
				}
			}
		}
	}
}

// Perform simple feature matching.  This just uses the SSD
// distance between two feature vectors, and matches a feature in the
// first image with the closest feature in the second image.  It can
// match multiple features in the first image to the same feature in
// the second image.
void ssdMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore) {
	int m = f1.size();
	int n = f2.size();

	matches.resize(m);
	totalScore = 0;

	double d;
	double dBest;
	int idBest;

	for (int i = 0; i<m; i++) {
		dBest = 1e100;
		idBest = 0;

		for (int j = 0; j<n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);

			if (d < dBest) {
				dBest = d;
				idBest = f2[j].id;
			}
		}	
		matches[i].id1 = f1[i].id;
		matches[i].id2 = idBest;
		matches[i].score = dBest;
		totalScore += matches[i].score;
	}
}

// TODO: Write this function to perform ratio feature matching.  
// This just uses the ratio of the SSD distance of the two best matches as the score
// and matches a feature in the first image with the closest feature in the second image.
// It can match multiple features in the first image to the same feature in
// the second image.  (See class notes for more information, and the sshMatchFeatures function above as a reference)
void ratioMatchFeatures(const FeatureSet &f1, const FeatureSet &f2, vector<FeatureMatch> &matches, double &totalScore)
{
	int m = f1.size();
	int n = f2.size();

	totalScore = 0;

	double d;
	double d_best;
	double d_2nd_best;

	int id_best;
	int id_2nd_best;

	for (int i = 0; i < m; i++) {
		d_best = 1e100;
		id_best = 0;
		d_2nd_best = 1e100;
		id_2nd_best = 0;
		for (int j = 0; j < n; j++) {
			d = distanceSSD(f1[i].data, f2[j].data);
			if (d < d_best) {
				d_2nd_best = d_best;
				id_2nd_best = id_best;
				d_best = d;
				id_best = f2[j].id;
			}
			else if (d < d_2nd_best && d >= d_best) {
				d_2nd_best = d;
				id_2nd_best = f2[j].id;
			}
		}
		FeatureMatch good_feature;
		good_feature.id1 = f1[i].id;
		good_feature.id2 = id_best;	
		good_feature.score = d_best / d_2nd_best;
		matches.push_back(good_feature);

		// cout << "match: " << matches[i].id1 << " " << matches[i].id2 << endl;
	}
}

// Convert Fl_Image to CFloatImage.
bool convertImage(const Fl_Image *image, CFloatImage &convertedImage) {
	if (image == NULL) {
		return false;
	}

	// Let's not handle indexed color images.
	if (image->count() != 1) {
		return false;
	}

	int w = image->w();
	int h = image->h();
	int d = image->d();

	// Get the image data.
	const char *const *data = image->data();

	int index = 0;

	for (int y = 0; y<h; y++) {
		for (int x = 0; x<w; x++) {
			if (d < 3) {
				// If there are fewer than 3 channels, just use the
				// first one for all colors.
				convertedImage.Pixel(x, y, 0) = ((uchar)data[0][index]) / 255.0f;
				convertedImage.Pixel(x, y, 1) = ((uchar)data[0][index]) / 255.0f;
				convertedImage.Pixel(x, y, 2) = ((uchar)data[0][index]) / 255.0f;
			}
			else {
				// Otherwise, use the first 3.
				convertedImage.Pixel(x, y, 0) = ((uchar)data[0][index]) / 255.0f;
				convertedImage.Pixel(x, y, 1) = ((uchar)data[0][index + 1]) / 255.0f;
				convertedImage.Pixel(x, y, 2) = ((uchar)data[0][index + 2]) / 255.0f;
			}

			index += d;
		}
	}

	return true;
}

// Convert CFloatImage to CByteImage.
void convertToByteImage(CFloatImage &floatImage, CByteImage &byteImage) {
	CShape sh = floatImage.Shape();

	assert(floatImage.Shape().nBands == byteImage.Shape().nBands);
	for (int y = 0; y<sh.height; y++) {
		for (int x = 0; x<sh.width; x++) {
			for (int c = 0; c<sh.nBands; c++) {
				float value = floor(255 * floatImage.Pixel(x, y, c) + 0.5f);

				if (value < byteImage.MinVal()) {
					value = byteImage.MinVal();
				}
				else if (value > byteImage.MaxVal()) {
					value = byteImage.MaxVal();
				}

				// We have to flip the image and reverse the color
				// channels to get it to come out right.  How silly!
				byteImage.Pixel(x, sh.height - y - 1, sh.nBands - c - 1) = (uchar)value;
			}
		}
	}
}

// Compute SSD distance between two vectors.
double distanceSSD(const vector<double> &v1, const vector<double> &v2) {
	int m = v1.size();
	int n = v2.size();

	if (m != n) {
		// Here's a big number.
		return 1e100;
	}

	double dist = 0;

	for (int i = 0; i<m; i++) {
		dist += pow(v1[i] - v2[i], 2);
	}


	return sqrt(dist);
}

// Transform point by homography.
void applyHomography(double x, double y, double &xNew, double &yNew, double h[9]) {
	double d = h[6] * x + h[7] * y + h[8];

	xNew = (h[0] * x + h[1] * y + h[2]) / d;
	yNew = (h[3] * x + h[4] * y + h[5]) / d;
}

// Compute AUC given a ROC curve
double computeAUC(vector<ROCPoint> &results)
{
	double auc = 0;
	double xdiff, ydiff;
	for (int i = 1; i < (int)results.size(); i++)
	{
		//fprintf(stream,"%lf\t%lf\t%lf\n",thresholdList[i],results[i].falseRate,results[i].trueRate);
		xdiff = (results[i].falseRate - results[i - 1].falseRate);
		ydiff = (results[i].trueRate - results[i - 1].trueRate);
		auc = auc + xdiff*results[i - 1].trueRate + xdiff*ydiff / 2;

	}
	return auc;
}