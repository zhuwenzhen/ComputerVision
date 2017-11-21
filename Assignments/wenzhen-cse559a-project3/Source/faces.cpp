
/////////////////////////////////////////////////////////////////////////////////////////////////
//	Project 4: Eigenfaces                                                                      //
//  CSE 455 Winter 2003                                                                        //
//	Copyright (c) 2003 University of Washington Department of Computer Science and Engineering //
//                                                                                             //
//  File: faces.cpp                                                                            //
//	Author: David Laurence Dewey                                                               //
//	Contact: ddewey@cs.washington.edu                                                          //
//           http://www.cs.washington.edu/homes/ddewey/                                        //
//                                                                                             //
/////////////////////////////////////////////////////////////////////////////////////////////////



#include "stdafx.h"

#include "jacob.h"

Faces::Faces()
:
Array<Face>(),
width(0),
height(0),
vector_size(0)
{
	//empty
}

Faces::Faces(int count, int width, int height)
:
Array<Face>(count),
width(width),
height(height),
vector_size(width*height)
{
	for (int i=0; i<getSize(); i++) {
		(*this)[i].resize(width, height, 1);
	}
}

void Faces::load(BinaryFileReader& file)
{
	resize(file.readInt());
	width=file.readInt();
	height=file.readInt();
	vector_size=width*height;
	for (int i=0; i<getSize(); i++) {
		(*this)[i].load(file);
	}
	average_face.load(file);
	std::cout << "Loaded faces from '" << file.getFilename() << "'" << std::endl;

}

void Faces::load(std::string filename)
{
	BinaryFileReader file(filename);
	load(file);
}

void Faces::save(std::string filename) const
{
	BinaryFileWriter file(filename);
	save(file);
}

void Faces::save(BinaryFileWriter& file) const
{
	file.write(getSize());
	file.write(width);
	file.write(height);
	for (int i=0; i<getSize(); i++) {
		(*this)[i].save(file);
	}
	average_face.save(file);
	std::cout << "Saved faces to '" << file.getFilename() << "'" << std::endl;
}

void Faces::output(std::string filepattern) const
{
	for (int i=0; i<getSize(); i++) {
		// normalize for output
		Image out_image;
		(*this)[i].normalize(0.0, 255.0, out_image);
		std::string filename=Functions::filenameNumber(filepattern, i, getSize()-1);
		out_image.saveTarga(filename);
	}
}

void Faces::eigenFaces(EigFaces& results, int n) const
{
	// size the results vector
	results.resize(n);
	results.setHeight(height);
	results.setWidth(width);

	// allocate matrices
	double **matrix = Jacobi::matrix(1, vector_size, 1, vector_size);
	double **eigmatrix = Jacobi::matrix(1, vector_size, 1, vector_size);
	double *eigenvec = Jacobi::vector(1, vector_size);
        
	// --------- TODO #1: fill in your code to prepare a matrix whose eigenvalues and eigenvectors are to be computed.
	// Also be sure you store the average face in results.average_face (A "set" method is provided for this).
	Face average;
	average.resize(width,height,1);

	for(int y=0;y<height;y++){
		for(int x=0;x<width;x++){
			average.pixel(x,y,0) = 0;
			for(int count=0;count<vector_size;count++){
				average.pixel(x,y,0) += (*this)[count].pixel(x,y,0);
			}
			average.pixel(x,y,0) = average.pixel(x,y,0)/vector_size;
		}
	}
	results.setAverage(average);

	//clear all the gradients
	for(int x=1;x<=vector_size;x++){
		for(int y=1;y<=vector_size;y++){
			matrix[x][y]=0;
		}
	}

	double **AT = Jacobi::matrix(1,vector_size,1,vector_size);
	double **A = Jacobi::matrix(1,vector_size,1,vector_size);
	double temp= 0 ;
	for (int count=1;count<=vector_size;count++){
		for(int y=0;y<height;y++){
			for(int x=0;x<width;x++){
				temp = (*this)[count-1].pixel(x,y,0) - average.pixel(x,y,0);
				A[x+y*width+1][count] = temp;
				AT[count][x+y*width+1] = temp;
			}
		}
	}
	for(int y=1;y<=vector_size;y++){
		for(int x=1;x<=vector_size;x++){
			//Jacobi matrix call element in different way, first for row num, second for col num, differ from pixel
			for(int ind=1;ind<=vector_size;ind++){
				matrix[y][x] += AT[y][ind] * A[ind][x] ;
			}
		}
	}

	// find eigenvectors
	int nrot;
	Jacobi::jacobi(matrix, vector_size, eigenvec, eigmatrix, &nrot);
	// sort eigenvectors
	Array<int> ordering;
	sortEigenvalues(eigenvec, ordering);
	for (int i=0; i<n; i++) {
		for (int k=0; k<vector_size; k++) {
			results[i][k] = eigmatrix[k+1][ordering[i]+1];
		}
	}

	double *length = Jacobi::vector(1, n);
	for (int i=0; i<n; i++) {
		for (int k=0; k<vector_size; k++) {
			length[i+1] += results[i][k] * results[i][k];
		}
	}
	for (int i=0; i<n ;i++){
		length[i+1] = sqrt (length[i+1]);
	}
	for (int i=0; i<n; i++) {
		for (int k=0; k<vector_size; k++) {
			results[i][k] = results[i][k] / length[i+1];
		}
	}

	// free matrices
	Jacobi::free_matrix(matrix, 1, vector_size, 1, vector_size);
	Jacobi::free_matrix(eigmatrix, 1, vector_size, 1, vector_size);
	Jacobi::free_vector(eigenvec, 1, vector_size);
}



int Faces::getWidth() const
{
	return width;
}

int Faces::getHeight() const
{
	return height;
}

void Faces::setWidth(int width)
{
	//width=width;
	vector_size=width*height;
}

void Faces::setHeight(int height)
{
	//height=height;
	vector_size=width*height;
}

void Faces::sortEigenvalues(double *eigenvec, Array<int>& ordering) const
{
	// for now use simple bubble sort
	ordering.resize(vector_size);
	std::list<EigenVectorIndex> list;
	int size = getSize();
	for (int i=0; i< size; i++) {
		EigenVectorIndex e;
		e.eigenvalue=eigenvec[i+1];
		e.index=i;
		list.push_back(e);
	}
	bool change=true;
	list.sort();
	std::list<EigenVectorIndex>::iterator it=list.begin();
	int n=0;
	while (it!=list.end()) {
		ordering[n] = (*it).index;
		it++;
		n++;
	}
}

