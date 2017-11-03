//============================================================================
// Name        : fun.cpp
// Author      : l
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <ctime>
#include <unistd.h>


#include "lib3ds/file.h"
#include "lib3ds/io.h"
#include "lib3ds/mesh.h"

#include "BVH.h"
#include "Triangle.h"
#include "ComplexVector3.h"

using std::vector;

using namespace std;

const double PI = 3.14159265358979323846;
const double c = 299792458;
const CXF j = CXF(0, 1);

//class RayIndex {
//public:
//	int ray_one;
//	int ray_two;
//	int ray_three;
//	int ray_four;
//	int ray_center;
//};
//struct Plane{
//    Plane() {}
//
//    Plane(const Vector3 & p, const Vector3 & n)
//        : p(p) // point
//        , n(n) //normal
//    {}
//
//    Vector3 p;
//    Vector3 n;
//};

string IntToStr(int n);

bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, double &t);
void ComplexField(const Vector3 &R, const double &kr, const double frequency, CXV3 &C);
double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2);
void rotate_THETA_then_PHI(Vector3* R1, int n, double THETA, double PHI);


double SafeAcos(double x);



void getMatlabData(const string &POL,  const int &PointOfAngle, double *RCSTheta, double *RCSPhi);


void initializeEF(const string &POL, const Vector3 &kin, Vector3 &theta_h, Vector3 &phi_h,
		Vector3 &Ei);


void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube,
		const Vector3 &tube0, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle);


void PO1CEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in,
		const Vector3 *k_tube,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle);





int GeometricalOptics(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ei,
		const int &MaximumNumberOfReflections,
		 Vector3 &kdir, Vector3 *PO,
		 double *kdotr, Vector3 *RefPoint);

void calAllGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
		const Vector3 &Ein,const int &MaximumNumberOfReflections,
		  Vector3 *k_tube,
		 int *hitnum_tube, Vector3 ** EpoA,
		double **kdotr, Vector3 **RefPoint);




void process_mem_usage(double& vm_usage, double& resident_set);
void getCSTformatData(const double *Theta, const double *Phi, const double *ReE_Theta, const double *ImE_Theta,
		const double *ReE_Phi, const double *ImE_Phi, const int &PointOfAngle, const int &NumberOfFrequencyPoints,
		const int ii);
void calEq24( Vector3 **EpoA,  double **kdotr,  Vector3 **RefPoint,
		const Vector3 &kin, const int *hitnum_tube, const int &MaximumNumberOfReflections,
		const int &tnum, const Vector3 &tube0, const double &cellSize, const double &k0,
		const Vector3 &theta_h, const Vector3 &phi_h, CXF *AthetaTotal,  CXF *AphiTotal,
		double** RetDis, CXV3**EcPO, CXF **Ex, CXF** Ey, CXF **SingleRayAtheta, CXF **SingleRayAphi);

void getSourcePoint(const double &cellSize, const int &tnum, const int xdnum, Vector3 *antt){
	Vector3 antc;

	for (int i = 0; i < tnum; i++) {
		antt[i].x = 0 + (i % xdnum) * cellSize;
		antt[i].y = 0 + i / xdnum * cellSize;
		antt[i].z = 0;
	}

//
//	for (int i = 0; i < rnum; i++) {
//		  cout<<antr[i].x<<","<<antr[i].y<<","<<antr[i].z<<"antr"<<endl;
//	}

	antc = (antt[0] + antt[tnum - 1]) / 2.0;



	for (int i = 0; i < tnum; i++) {
		antt[i] = antt[i] - antc;
	}
}
void movAntToFarField(const Vector3 *antt, const int &tnum,  double &distance, Vector3 *tube,  Vector3 &kin){
	Vector3  cros_1, cros_2, Vecotr3, ant_Nd;
	cros_1 = antt[1] - antt[0];
	cros_2 = antt[tnum-1] - antt[0];

	kin = normalize(cros_2 ^ cros_1);

	ant_Nd = distance * kin;


	for (int i = 0; i < tnum; i++) {

		tube[i] = antt[i] - ant_Nd;
	}

}

void rotateTHETAThenPHI(const Vector3 *R1, int n, double THETA, double PHI, Vector3 *R3) {
//	//Roate X Axis (THETA)
//	Vector3* R2 = new Vector3[n];
//	for (int i = 0; i < n; i++) {
//		R2[i].x = R1[i].x;
//		R2[i].y = cos(THETA) * R1[i].y - sin(THETA) * R1[i].z;
//		R2[i].z = sin(THETA) * R1[i].y + cos(THETA) * R1[i].z;
//	}
	//Roate Y Axis (THETA)
	Vector3* R2 = new Vector3[n];

	for (int i = 0; i < n; i++) {
		R2[i].x = cos(THETA) * R1[i].x - sin(THETA) * R1[i].z;
		R2[i].y = R1[i].y;
		R2[i].z = sin(THETA) * R1[i].x + cos(THETA) * R1[i].z;
	}

	//Roate Z Axis (PHI)
	for (int i = 0; i < n; i++) {
		R3[i].x = cos(PHI) * R2[i].x - sin(PHI) * R2[i].y;
		R3[i].y = sin(PHI) * R2[i].x + cos(PHI) * R2[i].y;
		R3[i].z = R2[i].z;
	}

	delete [] R2;
}


















void getAntPlane(const double&THETA, const double &PHI, const double &distance,
		const double &cellSize, const int &tnum, const int xdnum, Vector3 *antt,
		Vector3  &ant_N);

void getCSTformatDataf(const double *Theta, const double *Phi,  double **ReE_Theta,  double **ImE_Theta,
		 double **ReE_Phi, double  **ImE_Phi, const int &PointOfAngle, const int &NOFP){
    ofstream outFile;
    string CSTformat;




		cout<< CSTformat << "  \n";

		for(int m=0 ; m<NOFP ; m++){
				CSTformat="rcs" + IntToStr(m) +".dat";
				outFile.open(CSTformat.c_str());

			for(int n =0; n<PointOfAngle ; ++n){

				outFile << Theta[n]/PI*180<<"\t"<<Phi[n]/PI*180<<"\t"<<ReE_Theta[n][m]<<"\t"
					<< ImE_Theta[n][m]<<"\t"<<ReE_Phi[n][m]<<"\t"<<ImE_Phi[n][m]<<endl;


			}
	        outFile.close();



		}

}











int main(int argc, char *argv[]) {
	double START1, END1;


//-------------------------------input .3ds model---------------------------------
//	string input = "/home/user/cuda-workspace/lib3/Ball.3ds";
//	string input = "/home/user/cuda-workspace/lib3/plate1515.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 45d.3ds";
	string input = "/home/user/cuda-workspace/lib3/dihedral m.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 22pt5d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/trihedral30.3ds";
//	string input = "/home/user/cuda-workspace/lib3/threeinone.3ds";
	Lib3dsFile* fin = lib3ds_file_load(input.c_str());

	if (!fin) {
		cerr << "*****ERROR****\n Loading file failed" << input << endl;
		exit(EXIT_FAILURE);
	} else {
		cout << "Read ok" << endl;
	}
	Lib3dsMesh* mesh = fin->meshes;
	size_t num_face = mesh->faces;
//	size_t num_vert = mesh->points;

//	cout<<"#Face ="<<num_face<<endl;
//	cout<<"#Vertex="<<num_vert<<endl;

	Lib3dsPoint* vert = mesh->pointL;
	Lib3dsFace* face = mesh->faceL;

//	for(size_t i=0; i<num_vert; i++){
//		cout<<vert[i].pos[2]<<endl;
//	}
//	exit (0);

	vector<Object*> objects;

	for (size_t i = 0; i < num_face; i++) {
		objects.push_back(
				new Triangle(
						Vector3(vert[face[i].points[0]].pos[0],
								vert[face[i].points[0]].pos[1],
								vert[face[i].points[0]].pos[2]),
						Vector3(vert[face[i].points[1]].pos[0],
								vert[face[i].points[1]].pos[1],
								vert[face[i].points[1]].pos[2]),
						Vector3(vert[face[i].points[2]].pos[0],
								vert[face[i].points[2]].pos[1],
								vert[face[i].points[2]].pos[2])));
	}
	BVH bvh(&objects);
//-------------------------------input .3ds model---------------------------------








//   ------------------------------input Angle file-----------------------------------------
//	string ExternalAngleFile = "N";
//    std::vector<double> phi2;
//
//
//	if (ExternalAngleFile == "Y") {
//		std::ifstream ifile("/home/user/cuda-workspace/lib3/phi2.txt", std::ios::in);
//
//
//	    //check to see that the file was opened correctly:
//	    if (!ifile.is_open()) {
//	        std::cerr << "There was a problem opening the input file!\n";
//	        exit(1);//exit or do additional error checking
//	    }
//
//	    double num = 0.0;
//	    //keep storing values from the text file so long as data exists:
//	    while (ifile >> num) {
//	        phi2.push_back(num);
//	    }
//	    //verify that the scores were stored correctly:
//	//	    for (int i = 0; i < phi2.size(); ++i) {
//	//	        cout << phi2[i] << std::endl;
//	//	    }
//
//	} else {
////		double Pmin = -45;
////		double Pmax = 45;
////		double Pstep = 1;
////		int PointOfAngle = (Pmax-Pmin)/Pstep + 1;
//
//
//	}



//	std::ifstream ifile("/home/user/cuda-workspace/lib3/phi2.txt", std::ios::in);
//    std::vector<double> phi2;
//
//    //check to see that the file was opened correctly:
//    if (!ifile.is_open()) {
//        std::cerr << "There was a problem opening the input file!\n";
//        exit(1);//exit or do additional error checking
//    }
//
//    double num = 0.0;
//    //keep storing values from the text file so long as data exists:
//    while (ifile >> num) {
//        phi2.push_back(num);
//    }

    //verify that the scores were stored correctly:
//	    for (int i = 0; i < phi2.size(); ++i) {
//	        cout << phi2[i] << std::endl;
//	    }
//   ------------------------------input Angle file----------------------------------------




		double fcell = 15e9;					//max sweep frequecny
		double f ;

		double Fmin = 15e9;
		double Fmax = 15e9;
		double Fstep = 1e9;
		int NumberOfFrequencyPoints = (Fmax-Fmin)/Fstep+1;
		double lambdac = c / fcell;
//		cout<<"L="<<lambda<<endl;
//		exit(0);

		double Lx = 0.6;
		double Ly = 0.6;
		double rayPerWavelenth = 16;
		double cellSize = lambdac/rayPerWavelenth;
		double Pmin = -45;
		double Pmax = 45;
		double Pstep = 0.5;
		int PointOfAngle = (Pmax-Pmin)/Pstep + 1;
		int xdnum = Lx / cellSize;
		int ydnum = Ly / cellSize;
		const int mhpd = 3; //Maximum hit-point depth or MaximumNumberOfReflections

		double distance = 40;
//		distance > 2D^2/lambda
		string POL = "H";
//		int rnum = xdnum * ydnum;
//		int tnum = (xdnum - 1) * (ydnum - 1);

//		int rnum = (xdnum+1) * (ydnum+1);
		int tnum = xdnum * ydnum;
		Vector3* source = new Vector3[tnum]; //raytube
		Vector3* tube = new Vector3[tnum]; //raytube
		Vector3* k_tube = new Vector3[tnum];
		Vector3** EpoA = new Vector3*[tnum];
		  for(int i = 0; i < tnum; ++i){
			  EpoA[i] = new Vector3[mhpd];
		  }

		Vector3* kt_in = new Vector3[tnum];
		Vector3* refpoint_tube = new Vector3[tnum]; //The source point on the object
		Vector3** RefPoint = new Vector3*[tnum]; //The source point on the object
		  for(int i = 0; i < tnum; ++i){
			  RefPoint[i] = new Vector3[mhpd];
		  }

		int *hitnum_tube = new int[tnum];

		double** kdotr = new double*[tnum];  //Point of emission to the point of reflection distance
		  for(int i = 0; i < tnum; ++i){
			  kdotr[i] = new double[mhpd];
		  }



//		Vector3* retpoint_tube = new Vector3[tnum]; //The Return point on incident plane


		double* RCSTheta = new double[PointOfAngle];
		double* RCSPhi =new double[PointOfAngle];

		double* PHI = new  double[PointOfAngle];
		double* THETA = new double[PointOfAngle];

		double** ReE_Theta = new double*[PointOfAngle];
		double** ImE_Theta = new double*[PointOfAngle];
		double** ReE_Phi = new double*[PointOfAngle];
		double** ImE_Phi = new double*[PointOfAngle];
		 for(int i = 0; i < PointOfAngle; ++i){

			 ReE_Theta[i] = new double[NumberOfFrequencyPoints];
			 ImE_Theta[i] = new double[NumberOfFrequencyPoints];
			 ReE_Phi[i] = new double[NumberOfFrequencyPoints];
			 ImE_Phi[i] =new double[NumberOfFrequencyPoints];

		}

		memset(hitnum_tube, 0, tnum);
		memset(kt_in, 0, tnum);

		double** RetDis = new double*[tnum]; //Return Distance
		for(int i = 0; i < tnum; i++){
			RetDis[i] = new double[mhpd];
		}
		CXF* Atheta_total = new CXF[mhpd+1];
		CXF* Aphi_total = new CXF[mhpd+1];

		CXV3** EcPO = new CXV3*[tnum];
		for(int i = 0; i < tnum; i++){
			EcPO[i] = new CXV3[mhpd];
		 }
		CXF** Ex = new CXF*[tnum];
		CXF** Ey = new CXF*[tnum];
		CXF** SingleRayAtheta  = new CXF*[tnum];
		CXF** SingleRayAphi =new CXF*[tnum];
		for(int i = 0 ; i<tnum ; i++){
			Ex[i] = new CXF[mhpd];
			Ey[i] = new CXF[mhpd];
			SingleRayAtheta[i] = new CXF[mhpd];
			SingleRayAphi[i] =  new CXF[mhpd];

		}

		START1=clock();


//		double PHI;
		Vector3 phi_h;
		Vector3 theta_h;
		Vector3 Ei;

		getSourcePoint(cellSize, tnum, xdnum, source);

//		for (int i = 0  ; i < phi2.size(); i++) {
		for (int i = 0  ; i < PointOfAngle; i++) {

			PHI[i] = (i*Pstep + Pmin) * PI / 180;
//			PHI[i] = phi2[i] * PI / 180;

			THETA[i] = -90 * PI / 180;


	//		PHI =  phi2[i]*PI/180;
			cout<<"PHI = "<<PHI[i]/PI*180<<endl;
//			cout<<"PHI = "<<PHI[i]<<endl;
			rotateTHETAThenPHI(source, tnum, THETA[i], PHI[i], tube) ;

			movAntToFarField(tube, tnum, distance, tube, kt_in[0]);
//			getAntPlane(THETA[i], PHI[i], distance,
//					cellSize, tnum, xdnum, tube,
//					  kt_in[0]);



			initializeEF(POL, kt_in[0], theta_h, phi_h, Ei);



	//		cout<< "THETA_"<<"("<<theta_h.x<<","<<theta_h.y<<","<<theta_h.z<<")"<<endl;
	//		cout<< "PHI_"<<"("<<phi_h.x<<","<<phi_h.y<<","<<phi_h.z<<")"<<endl;
			calAllGeoOpt(tnum, tube ,bvh, kt_in[0], Ei, mhpd,
						  k_tube,  hitnum_tube, EpoA,
						 kdotr, RefPoint);

//			exit(0);
			for(int ii = 0 ;ii < NumberOfFrequencyPoints; ii++){
				f = (ii*Fstep + Fmin);
				double k0 = 2 * PI * f * pow(c, -1);
				calEq24(EpoA, kdotr, RefPoint,
						kt_in[0], hitnum_tube, mhpd,
						tnum, tube[0], cellSize, k0,
						theta_h, phi_h, Atheta_total,  Aphi_total,
						 RetDis, EcPO, Ex, Ey, SingleRayAtheta,  SingleRayAphi);




	//			exit(0);

				Atheta_total[mhpd]=CXF(0,0);
				Aphi_total[mhpd]=CXF(0,0);





				for(int n =0; n< mhpd; n++){

					Atheta_total[mhpd]+=Atheta_total[n];
					Aphi_total[mhpd]+=Aphi_total[n];

				}



			 ReE_Theta[i][ii] = Atheta_total[mhpd].real();
			 ReE_Phi[i][ii] = Aphi_total[mhpd].real();
			 ImE_Phi[i][ii] = Aphi_total[mhpd].imag();
			 ImE_Theta[i][ii] = Atheta_total[mhpd].imag();

	//
				RCSTheta[i] = 10* log10(PI * 4 * abs(Atheta_total[mhpd]) * abs(Atheta_total[mhpd]));
				RCSPhi[i]= 10 * log10(PI * 4 * abs(Aphi_total[mhpd]) * abs(Aphi_total[mhpd]));




	//			RCSABS[i] = 10 * log10(PI * 4 * abs(Aphi_total[3]) * abs(Aphi_total[3])+
	//					PI * 4 * abs(Atheta_total[3]) * abs(Atheta_total[3]));
	//			RCSABS[i] =Aphi_total[1];
		//
				cout << "RCS(theta) ="
						<< RCSTheta[i]
						<< " dBsm" << endl;
				cout << "RCS(phi) ="<<
						RCSPhi[i]<< " dBsm"
						<< endl;

	//			exit(0);

//				getCSTformatData(THETA, PHI, ReE_Theta, ImE_Theta,
//								ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints, ii);
//				getCSTformatDataf(THETA, PHI, ReE_Theta, ImE_Theta,
//								ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints, ii, f);




	//			getAngleData(PointOfAngle, RCSABS);

			}
//			getCSTformatDataf(THETA, PHI, ReE_Theta, ImE_Theta,
//							ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints, f);
			getCSTformatDataf(THETA, PHI, ReE_Theta, ImE_Theta,
							ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints);
//			END1=clock();
//				cout << endl << "進行運算所花費的時間：" << (END1- START1) / CLOCKS_PER_SEC << " S" << endl;
//				exit(0);
//


		}
//		getCSTformatData(THETA, PHI, ReE_Theta, ImE_Theta,
//						ReE_Phi, ImE_Phi, phi2.size(), NumberOfFrequencyPoints, ii);
//		getCSTformatDataf(THETA, PHI, ReE_Theta, ImE_Theta,
//						ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints);




//		getMatlabData(POL, PointOfAngle, RCSTheta, RCSPhi);
		delete[] source;

		delete[] tube;
		delete[] k_tube;
		delete[] hitnum_tube;

		for(int i = 0; i<tnum ; ++i){

		    delete [] EpoA[i];
		  }
		delete[] EpoA;

		for(int i = 0; i < tnum; i++){
			delete RetDis[i];
			delete Ex[i];
			delete Ey[i];
			delete SingleRayAtheta[i];
			delete SingleRayAphi[i];
			delete EcPO[i];

		}

		delete[] RetDis;
		delete[] EcPO;

		delete []Ex;
		delete []Ey;
		delete []SingleRayAtheta;
		delete []SingleRayAphi;


		for(int i = 0; i<tnum ; ++i){

		    delete [] kdotr[i];
		  }
		delete[] kdotr;


		delete[] refpoint_tube; //The source point on the object

		for(int i = 0; i<tnum ; ++i){

		    delete [] RefPoint[i];
		  }
		delete [] RefPoint;

//		delete[] Atheta;
//		delete[] Aphi;
		delete[] kt_in;
		delete[] Atheta_total;
		delete[] Aphi_total;
		delete[] RCSTheta;
		delete[] RCSPhi;
		delete[] PHI;
		delete[] THETA;
		 for(int i = 0; i <PointOfAngle ; ++i){

			delete[] ReE_Theta[i];
			delete[] ImE_Theta[i];
			delete[] ReE_Phi[i];
			delete[] ImE_Phi[i];

		}

		delete[] ReE_Theta;
		delete[] ImE_Theta;
		delete[] ReE_Phi;
        delete[] ImE_Phi;


	cout << (double)clock() / CLOCKS_PER_SEC << " S"<<"\n";
	   double vm, rss;

	   process_mem_usage(vm, rss);
	   cout << "VM: " << vm << "; RSS: " << rss << endl;

//	}   frequency loop
//    ofstream outFile;
//    int Number_of_files=2;
//    string CSTformat;


//  for (int iii=0;iii<Number_of_files;iii++)
//  {
//        CSTformat="rcs" + IntToStr(iii) +".dat";
//        cout<< CSTformat << "  \n";
//
//        outFile.open(CSTformat.c_str());
//        outFile << CSTformat<<" : Writing this to a file.\n\n"
//        		"Fake CST Farfield Format V1 \n"
//
//        		"Dimension     = 1\n"
//        		"Frequency     = 1.5e+010\n"
//        		"Samples       = 181\n"
//        		"WaveAmplitude = 1\n"
//        		"Type          = MONOSTATICRCS\n\n"
//
//        		"// = Theta Phi Re(E_Theta) Im(E_Theta) Re(E_Phi) Im(E_Phi)\n";
//
//        outFile.close();
//  }
	return 0;

}






//Ray-Plane Intersection
//Ref: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-plane-and-ray-disk-intersection
bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, double &t) {
	// assuming vectors are all normalized
	double denom = n * l;
	if (denom > 1e-15) {
		Vector3 p0l0 = p0 - l0;
		t = p0l0 * n / denom;
		return (t >= 0);
	}

	return false;
}
void ComplexField(const Vector3 &R, const double &kr, const double k0, CXV3 &C) {
//	CXF expvalue = exp(j * (kr) / c * frequency * multiply2 * PI);

	CXF expvalue = exp(j * (kr) *k0);
	C = CXV3(CXF(R.x, 0) * expvalue, CXF(R.y, 0) * expvalue,
			CXF(R.z, 0) * expvalue);
}

double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2) {
	double l1 = length(v1);
	double l2 = length(v2);
	double angle = 0;
//	cout << v2.x << "," << v2.y << "," << v2.z << endl;


	angle = SafeAcos(v1 * v2 / (l1 * l2));
	angle = angle * 180 / PI;
//	cout << angle << endl;

	return angle;
}

double SafeAcos(double x) {
	if (x < -1.0)
		x = -1.0;
	else if (x > 1.0)
		x = 1.0;
	return acos(x);
}

void rotate_THETA_then_PHI(Vector3 *R1, int n, double THETA, double PHI) {
//	//Roate X Axis (THETA)
//	Vector3* R2 = new Vector3[n];
//	for (int i = 0; i < n; i++) {
//		R2[i].x = R1[i].x;
//		R2[i].y = cos(THETA) * R1[i].y - sin(THETA) * R1[i].z;
//		R2[i].z = sin(THETA) * R1[i].y + cos(THETA) * R1[i].z;
//	}
	//Roate Y Axis (THETA)
	Vector3* R2 = new Vector3[n];

	for (int i = 0; i < n; i++) {
		R2[i].x = cos(THETA) * R1[i].x - sin(THETA) * R1[i].z;
		R2[i].y = R1[i].y;
		R2[i].z = sin(THETA) * R1[i].x + cos(THETA) * R1[i].z;
	}

	//Roate Z Axis (PHI)
	for (int i = 0; i < n; i++) {
		R1[i].x = cos(PHI) * R2[i].x - sin(PHI) * R2[i].y;
		R1[i].y = sin(PHI) * R2[i].x + cos(PHI) * R2[i].y;
		R1[i].z = R2[i].z;
	}

	delete [] R2;
}




void getAngleData(const int &PointOfAngle, const CXF *OneAngleNumberOfHit){
	ofstream angledata("angledata.out");
	for (int i = 0; i < PointOfAngle; i++) {
		 angledata <<  OneAngleNumberOfHit[i] <<endl;

		ifstream ifile("angledata.out");
	}
	angledata.close();
}

void getMatlabData(const string &POL,  const int &PointOfAngle, double *RCSTheta, double *RCSPhi){

		if(POL =="V"){

		ofstream file("ThetapolRcsTheta.out");
		for (int i = 0; i < PointOfAngle; i++) {
			file << RCSTheta[i] <<endl;
			ifstream ifile("ThetapolRcsTheta.out");
		}
		file.close();

		ofstream file1("ThetapolRcsPhi.out");
		for (int i = 0; i < PointOfAngle; i++) {
			file1 << RCSPhi[i] << endl;
			ifstream ifile1("ThetaPolRcsPhi.out");
		}
		file1.close();
		}else{
		ofstream file("PhipolRcsTheta.out");
		for (int i = 0; i < PointOfAngle; i++) {
			file << RCSTheta[i] <<endl;
			ifstream ifile("PhipolRcsTheta.out");
		}
		file.close();

		ofstream file1("PhipolRcsPhi.out");
		for (int i = 0; i < PointOfAngle ; i++) {
			file1 << RCSPhi[i] << endl;
			ifstream ifile1("PhipolRcsPhi.out");
		}
		file1.close();
		}
}



void initializeEF(const string &POL, const Vector3 &kin, Vector3 &theta_h, Vector3 &phi_h,
		Vector3 &Ei){

	Vector3 Nz = Vector3(0, 0, 1);

	 phi_h = kin ^ Nz;
	 theta_h = kin ^ phi_h;
//		cout << "phi_h = " << phi_h.x << "," << phi_h.y << "," << phi_h.z
//				<< endl;
//		cout << "theta_h = " << theta_h.x << "," << theta_h.y << ","
//				<< theta_h.z << endl;

	if (POL == "V") {
		Ei = theta_h; //v ->Y'
	} else {
		Ei = phi_h;   //h ->X'
	}
	//	 Ei = (Ei * theta_h) * theta_h + (Ei * phi_h) * phi_h;
}

void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube,
		const Vector3 &tube0, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle){


	for (int n = 0; n < tnum; n++) {

		if(hitnum_tube[n]!=0){
			IR_Angle[n] = Vector_Angle_d(kt_in[0], k_tube[n]);
//				cout<<IR_Angle[n]<<endl;

		}else{
			IR_Angle[n] = -1;
		}
		if (IR_Angle[n] >= 0 ){
			intersectPlane((-1 * kt_in[0]), tube0,
					refpoint_tube[n], (-1 * kt_in[0]), kdr_tuber[n]);


//			retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
			total_path_tube[n] = kdr_tuber[n] + kdr_tube[n];
//			total_path_tube[n] = 2*kdr_tube[n];
//			cout<<kdr_tube[n]<<"kdr_tube"<<endl;
//			cout<< kdr_tuber[n]<<"kdr_tbr"<<endl;

//			cout << total_path_tube[n] << " totalpath" << endl;

			ComplexField(E_tube[n], total_path_tube[n], f, Ec_tube[n]);
		} else {
			Ec_tube[n] = CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));
		}
		//		cout << "rethit = " << return_hit[n] << endl;
	}

}
void PO1CEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in,
		const Vector3 *k_tube,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle){

//	int a, b;

	for (int n = 0; n < tnum; n++) {

		if(hitnum_tube[n]!=0){
			IR_Angle[n] = Vector_Angle_d(kt_in[0], k_tube[n]);
//				cout<<IR_Angle[n]<<endl;

		}else{
			IR_Angle[n] = -1;
		}
		if (IR_Angle[n] > 0 ){
//			intersectPlane((-1 * kt_in[0]), tube[0],
//					refpoint_tube[n], (-1 * kt_in[0]), kdr_tuber[n]);
//			intersectPlane((-1 * kt_in[0]), tube[0],
//					refpoint_tube[n], k_tube[n], kdr_tuber[n]);

//			retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
//			total_path_tube[n] = kdr_tuber[n] + kdr_tube[n];
			total_path_tube[n] = 2*kdr_tube[n];
//			cout<<kdr_tube[n]<<"kdr_tube"<<endl;
//			cout<< kdr_tuber[n]<<"kdr_tbr"<<endl;

//			cout << total_path_tube[n] << " totalpath" << endl;

			ComplexField(E_tube[n], total_path_tube[n], f, Ec_tube[n]);
		} else {
			Ec_tube[n] = CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));
		}
		//		cout << "rethit = " << return_hit[n] << endl;
	}

}

void getAntPlane(const double&THETA, const double &PHI, const double &distance,
		const double &cellSize, const int &tnum, const int xdnum, Vector3 *antt,
		Vector3  &ant_N){
	Vector3 antc, cros_1, cros_2, ant_Nd;


	for (int i = 0; i < tnum; i++) {
		antt[i].x = 0 + (i % xdnum) * cellSize;
		antt[i].y = 0 + i / xdnum * cellSize;
		antt[i].z = 0;
	}

//
//	for (int i = 0; i < rnum; i++) {
//		  cout<<antr[i].x<<","<<antr[i].y<<","<<antr[i].z<<"antr"<<endl;
//	}

	antc = (antt[0] + antt[tnum - 1]) / 2.0;



	for (int i = 0; i < tnum; i++) {
		antt[i] = antt[i] - antc;
	}

	rotate_THETA_then_PHI(antt, tnum, THETA, PHI);


	cros_1 = antt[1] - antt[0];
	cros_2 = antt[tnum-1] - antt[0];

	ant_N = cros_2 ^ cros_1;
	ant_N = normalize(ant_N);
	ant_Nd = distance * ant_N;


	for (int i = 0; i < tnum; i++) {

		antt[i] = antt[i] - ant_Nd;
	}


//	antc = antc - ant_Nd;

}





void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;

}

int GeometricalOptics(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ei,
		const int &MaximumNumberOfReflections, Vector3 &kdir, Vector3 *PO,
		  double *kdotr, Vector3 *RefPoint) {

	Vector3 phi_h, theta_h, Xc, Yc, Zc, theta_ci_h, phi_ci_h, theta_cr_h,
			phi_cr_h, Nz, Ed, m_h, Er, nor, kipo;
	double theta_ci;
//	float phi_ci;
	int hitnumber = 0;            //
    double t=0;
	Ray ray(raypoint, kin);
	IntersectionInfo I;
//	for(int i = 0 ; i <MaximumNumberOfReflections; i++){
//		PO[i]=0*Ei;
//		kdotr[i] =0;
//		RefPoint[i] =Vector3(0,0,0);
//	}


	bool hit = bvh.getIntersection(ray, &I, false);

	while (hit) {
		hitnumber = hitnumber + 1;
//		cout<<hitnumber<<endl;
//		cout << "hitpoint = " << I.hit.x << "," << I.hit.y << "," << I.hit.z<< endl;
//      cout<<I.t<<"I.t"<<endl;
		nor = I.object->getNormal(I);
//      cout<<nor.x<<","<<nor.y<<","<<nor.z<<"nor"<<endl;

		//R1 formula (11)
		theta_ci = acos((-1 * ray.d) * nor);
//		phi_ci = 0;
//		cout<<"theta_ci = "<<theta_ci<<endl;

		if (theta_ci > 1e-5) {
			//R1 formula (6)
			m_h = (-1 * ray.d) ^ nor / sin(theta_ci);
//		    cout<<m_h.x<<","<<m_h.y<<","<<m_h.z<<"m"<<endl;

			//R1 formula (10)
			Xc = m_h ^ nor;
			Zc = -1 * nor;
			Yc = -1 * m_h;

			//R1 formula (12)
			theta_ci_h = Xc * cos(theta_ci) - Zc * sin(theta_ci);
			phi_ci_h = Yc;
			//R1 formula (13)
			theta_cr_h = (-1 * Xc * cos(theta_ci)) - Zc * sin(theta_ci);
			phi_cr_h = Yc;

			// get a reflection vector
			//using https://math.stackexchange.com/questions/13261/how-to-get-a-reflection-vector
			kipo = ray.d;
			kdir = ray.d - 2 * (ray.d * nor) * nor;

			//R1 formula (8) and (9)
			if (hitnumber == 1) {
//				Ed = (Ei * phi_ci_h) * phi_ci_h
//						+ (Ei * theta_ci_h) * theta_ci_h;
//				Ed = Ei;
//				Er = (-1) * (Ei * phi_ci_h) * phi_cr_h
//						+ (-1) * (Ei * theta_ci_h) * theta_cr_h;
				Er = (-1) * (Ei * phi_ci_h) * phi_cr_h
						- (-1) *(Ei * theta_ci_h) * theta_cr_h;
//				Er = (-1) * (Ei)+(2*(nor^Ei)^nor);


				PO[hitnumber-1]  = nor^(kipo^Ei)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t;
				RefPoint[hitnumber-1] = I.hit;


			} else {

				Ed = Er;
//				Er = (-1) * (Ed * phi_ci_h) * phi_cr_h
//						+ (-1) * (Ed * theta_ci_h) * theta_cr_h;
				Er = (-1) * (Ed * phi_ci_h) * phi_cr_h
						- (-1) * (Ed * theta_ci_h) * theta_cr_h;
//				Er = (-1) * (Ed)+(2*(nor^Ed)^nor);

				PO[hitnumber-1]  = nor^(kipo^Ed)/abs((kipo)*nor);

				kdotr[hitnumber-1] = I.t+t;
				RefPoint[hitnumber-1] = I.hit;

			}

		} else {
			kipo = ray.d;
			kdir = (-1) * ray.d;
			if (hitnumber == 1) {
//				Ed = Ei;
				Er =(-1)* Ei;
				PO[hitnumber-1] = nor^(kipo^Ei)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t;
				RefPoint[hitnumber-1] = I.hit;
			} else{
				Ed = Er;
				Er = (-1)*Ed;
				PO[hitnumber-1]  = nor^(kipo^Ed)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t + t;
				RefPoint[hitnumber-1] = I.hit;

			}
		}
		t += I.t;

		ray = Ray(I.hit, kdir);
		hit = bvh.getIntersection(ray, &I, false);
		if ((!hit) || (hitnumber == MaximumNumberOfReflections)) {
			break;

		}

	}
//		Eout = Er;

	return hitnumber;
}



void calAllGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
		const Vector3 &Ein,const int &MaximumNumberOfReflections,
		  Vector3 *k_tube, int *hitnum_tube, Vector3 ** EpoA,
		 double **kdotr, Vector3 **RefPoint){
	for (int n = 0; n < tnum; n++) {
		hitnum_tube[n] = GeometricalOptics(tube[n], bvh,  kin, Ein, MaximumNumberOfReflections,
				 k_tube[n], EpoA[n],
				 kdotr[n], RefPoint[n]);

	}

}


void getCSTformatData(const double *Theta, const double *Phi, const double *ReE_Theta, const double *ImE_Theta,
		const double *ReE_Phi, const double *ImE_Phi, const int &PointOfAngle, const int &NumberOfFrequencyPoints,
		const int ii){
	    ofstream outFile;
	    string CSTformat;

		CSTformat="rcs" + IntToStr(ii) +".dat";
		cout<< CSTformat << "  \n";
		outFile.open(CSTformat.c_str());

		for(int n=0 ; n<PointOfAngle ; n++){

			outFile << Theta[n]/PI*180<<"\t"<<Phi[n]/PI*180<<"\t"<<ReE_Theta[n]<<"\t"
					<< ImE_Theta[n]<<"\t"<<ReE_Phi[n]<<"\t"<<ImE_Phi[n]<<endl;


		}


	        outFile.close();

}
string IntToStr(int n)
{
    stringstream result;
    result << n;
    return result.str();
}

void calEq24( Vector3 **EpoA,  double **kdotr,  Vector3 **RefPoint,
		const Vector3 &kin, const int *hitnum_tube, const int &MaximumNumberOfReflections,
		const int &tnum, const Vector3 &tube0, const double &cellSize, const double &k0,
		const Vector3 &theta_h, const Vector3 &phi_h, CXF *AthetaTotal,  CXF *AphiTotal,
		double** RetDis, CXV3**EcPO, CXF **Ex, CXF** Ey, CXF **SingleRayAtheta, CXF **SingleRayAphi){

	for(int i = 0 ; i<MaximumNumberOfReflections ; i++){
		AthetaTotal[i]=CXF(0,0);
		AphiTotal[i] =CXF(0,0);
		for(int n =0 ; n<tnum ; n++){
			if(hitnum_tube[n]!=0 && hitnum_tube[n]>= i+1){
				intersectPlane((-1 * kin), tube0,
						RefPoint[n][i], (-1 * kin), RetDis[n][i]);

				RetDis[n][i]= RetDis[n][i] + kdotr[n][i];

				ComplexField(EpoA[n][i], RetDis[n][i], k0, EcPO[n][i]);


			}else{
				EcPO[n][i]=CXV3((CXF(0, 0)), CXF(0, 0), CXF(0, 0));

			}
			Ex[n][i] = EcPO[n][i] * phi_h;
			Ey[n][i] = EcPO[n][i] * theta_h;
			SingleRayAtheta[n][i] = j*k0/(2*PI)*( Ey[n][i] )  * cellSize * cellSize;
			SingleRayAphi[n][i] = j*k0/(2*PI)*(Ex[n][i] )  * cellSize * cellSize;
			AthetaTotal[i] += SingleRayAtheta[n][i];

			AphiTotal[i] += SingleRayAphi[n][i];


		}
	}

}




