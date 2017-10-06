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

class RayIndex {
public:
	int ray_one;
	int ray_two;
	int ray_three;
	int ray_four;
	int ray_center;

};
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
string IntToStr(int n)
{
    stringstream result;
    result << n;
    return result.str();
}

bool intersectPlane(const Vector3 &n, const Vector3 &p0, const Vector3 &l0,
		const Vector3 &l, double &t);
void ComplexField(const Vector3 &R, const double &kr, const double frequency, CXV3 &C);
double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2);
void rotate_THETA_then_PHI(Vector3* R1, int n, double THETA, double PHI);

void genrateAperture(const double&THETA, const double &PHI, const double &distance,
		const double cellSize, const int rnum,const int tnum,
		const double &Lx, const double &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
		Vector3  &ant_N);

double SafeAcos(double x);

void calEquationTwentySeven(const CXV3 *E, const double &k0,
		const double &cellSize, const int &tnum, const Vector3 &theta_h,
		const Vector3 &phi_h, const double &PHI, const double *IR_Angle, Vector3 *retpoint_antr,
		CXF &Atheta, CXF &Aphi);

void getMatlabData(const string &POL,  const int &PointOfAngle, double *RCSTheta, double *RCSPhi);


void initializeEF(const string &POL, const Vector3 &kin, Vector3 &theta_h, Vector3 &phi_h,
		Vector3 &Ei);


void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle);

//void Raytube_Numbering(const int &tnum, const int &xdnum, const Vector3 *antr, const Vector3 *antt, RayIndex *index) ;

void getAngleData(const int &PointOfAngle, const CXF *OneAngleNumberOfHit);


//int GeometricalOptics(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ein,
//		const int &MaximumNumberOfReflections,
//		Vector3 &hitpoint, Vector3 &kdir, Vector3 &Eout, double &t, Vector3 &NormalCrossEpo2,
//		Vector3 &PO1, double &tPO1);

//void AllRayGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
//		const Vector3 &Ein,const int &MaximumNumberOfReflections,
//		 Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
//		double *kdr_tube, int *hitnum_tube, Vector3 *NormalCrossEpo2, Vector3 *Epo1,
//		double *kdrPO1);
void PO1CEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle);



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



int GeometricalOptics(const Vector3 &raypoint, const BVH &bvh, const Vector3 &kin, const Vector3 &Ein,
		const int &MaximumNumberOfReflections,
		Vector3 &hitpoint, Vector3 &kdir, Vector3 &Eout, double &t, Vector3 *PO,
		Vector3 &PO1, double &tPO1, double *kdotr, Vector3 *RefPoint) {

	Vector3 phi_h, theta_h, Ei, Xc, Yc, Zc, theta_ci_h, phi_ci_h, theta_cr_h,
			phi_cr_h, Nz, Ed, m_h, Er, nor, kipo;
	double theta_ci;
//	float phi_ci;
	int hitnumber = 0;            //
	Ei = Ein;
    t=0;
	Ray ray(raypoint, kin);
	IntersectionInfo I;
	for(int i = 0 ; i <MaximumNumberOfReflections; i++){
		PO[i]=0*Ein;
		kdotr[i] =0;
		RefPoint[i] =Vector3(0,0,0);
	}


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
				Ed = Ei;
//				Er = (-1) * (Ei * phi_ci_h) * phi_cr_h
//						+ (-1) * (Ei * theta_ci_h) * theta_cr_h;
				Er = (-1) * (Ei * phi_ci_h) * phi_cr_h
						- (-1) *(Ei * theta_ci_h) * theta_cr_h;
//				Er = (-1) * (Ei)+(2*(nor^Ei)^nor);

				PO1 = nor^(kipo^Ed)/abs((kipo)*nor);
				tPO1 = I.t;
				PO[hitnumber-1] =PO1 = nor^(kipo^Ed)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t;
				RefPoint[hitnumber-1] = I.hit;


			} else {

				Ed = Er;
//				Er = (-1) * (Ed * phi_ci_h) * phi_cr_h
//						+ (-1) * (Ed * theta_ci_h) * theta_cr_h;
				Er = (-1) * (Ed * phi_ci_h) * phi_cr_h
						- (-1) * (Ed * theta_ci_h) * theta_cr_h;
//				Er = (-1) * (Ed)+(2*(nor^Ed)^nor);

				PO[hitnumber-1] =PO1 = nor^(kipo^Ed)/abs((kipo)*nor);

				kdotr[hitnumber-1] = I.t+t;
				RefPoint[hitnumber-1] = I.hit;

			}

		} else {
			kipo = ray.d;
			kdir = (-1) * ray.d;
			if (hitnumber == 1) {
				Ed = Ei;
				Er =(-1)* Ed;
				PO1 = nor^(kipo^Ed)/abs(kipo*nor);
				tPO1 = I.t;
				PO[hitnumber-1] =PO1 = nor^(kipo^Ed)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t;
				RefPoint[hitnumber-1] = I.hit;


			} else{
				Ed = Er;
				Er = (-1)*Ed;
				PO[hitnumber-1] =PO1 = nor^(kipo^Ed)/abs((kipo)*nor);
				kdotr[hitnumber-1] = I.t + t;
				RefPoint[hitnumber-1] = I.hit;

			}
		}
		t += I.t;
//		cout<<"Einput = "<<Ein.x<<","<<Ein.y<<","<<Ein.z<<endl;
//		cout<<"PO1 = "<<PO[0].x<<","<<PO[0].y<<","<<PO[0].z<<endl;


//		cout<<"kinput = "<<kin.x<<","<<kin.y<<","<<kin.z<<endl;
//		cout<<"kouput = "<<kdir.x<<","<<kdir.y<<","<<kdir.z<<endl;
//		cout<<"ki = "<<kipo.x<<","<<kipo.y<<","<<kipo.z<<endl;
//		cout<<"n = "<<nor.x<<","<<nor.y<<","<<nor.z<<endl;

//		cout<<"Eipo = "<<Ed.x<<","<<Ed.y<<","<<Ed.z<<endl;


		ray = Ray(I.hit, kdir);
		hit = bvh.getIntersection(ray, &I, false);
		if ((!hit) || (hitnumber == MaximumNumberOfReflections)) {
			break;

		}

	}
//
//		if (hitnumber==1){
//         NormalCrossEpo2 = nor^(kdir^Er);
//         NormalCrossEpo2 = nor^((kin)^Ei)/((kin)*nor);
//	        PO1 = 0* Ei ;
//         NormalCrossEpo2 = nor^(Er^nor);

//         cout<<kdir*nor<<endl;
//			cout<<"PO1 = "<<PO[0].x<<","<<PO[0].y<<","<<PO[0].z<<endl;


//		cout<<"nx(Erxn) = "<<.x<<","<<NormalCrossEpo2.y<<","<<NormalCrossEpo2.z<<endl;

//		}
		if (hitnumber==0){

//         NormalCrossEpo2 = nor^(kipo^Ed);
//         NormalCrossEpo2 = nor^(kdir^Er);
		   PO1 = 0* Ei ;

//         NormalCrossEpo2 = nor^(Er^nor);

//         cout<<kdir*nor<<endl;
//			cout<<"PO1 = "<<PO[0].x<<","<<PO[0].y<<","<<PO[0].z<<endl;

//		cout<<"nx(Erxn) = "<<NormalCrossEpo2.x<<","<<NormalCrossEpo2.y<<","<<NormalCrossEpo2.z<<endl;

		}
//		if(hitnumber==3){
//					cout<<"Einput = "<<Ein.x<<","<<Ein.y<<","<<Ein.z<<endl;
//
//					cout<<"PO3 = "<<PO[2].x<<","<<PO[2].y<<","<<PO[2].z<<endl;
//
//		}

		if (hitnumber==2){
//				NormalCrossEpo2 = nor^(kdir^Er)/abs(kdir*nor);
//				NormalCrossEpo2 = nor^(kipo^Ed);
//				cout<<"tpo1"<<tPO1<<endl;



//         NormalCrossEpo2 = nor^((kdir)^Er)/((kdir)*nor);
//         NormalCrossEpo2 = nor^(Er^nor);

//         cout<<kdir*nor<<endl;
//		cout<<"PO2 = "<<PO[1].x<<","<<PO[1].y<<","<<PO[1].z<<endl;
//
////		cout<<"nx(EIxn) = "<<NormalCrossEpo2.x<<","<<NormalCrossEpo2.y<<","<<NormalCrossEpo2.z<<endl;
//
		}

//        cout<<"Er"<<Er.x<<","<<Er.y<<","<<Er.z<<endl;
//        cout<<"epo"<<NormalCrossEpo2.x<<","<<NormalCrossEpo2.y<<","<<NormalCrossEpo2.z<<endl;


		Eout = Er;
		hitpoint = I.hit;

//		cout<<I.t<<endl;
//		cout<<t<<endl;
//		cout<<Ei.x<<","<<Ei.y<<","<<Ei.z<<" = ei"<<endl;
//		cout << "hitpoint = " << hitpoint.x << "," << hitpoint.y << ","
//				<< hitpoint.z << endl;

//		    cout<<Ed.x<<","<<Ed.y<<","<<Ed.z<<"ed"<<endl;
//		cout << Er.x << "," << Er.y << "," << Er.z << " =er " << endl;
//		cout << NormalCrossEpo2.x << "," << NormalCrossEpo2.y << "," << NormalCrossEpo2.z << " =enor " << endl;

//		cout << hitnumber << "hitnumber" << endl;

	return hitnumber;
}

//void AllRayGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
//		const Vector3 &Ein,const int &MaximumNumberOfReflections,
//		 Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
//		double *kdr_tube, int *hitnum_tube, Vector3 *NormalCrossEpo2, Vector3 *Epo1,
//		double *tPO1){
//	for (int n = 0; n < tnum; n++) {
//		hitnum_tube[n] = GeometricalOptics(tube[n], bvh,  kin, Ein, MaximumNumberOfReflections,
//				refpoint_tube[n], k_tube[n], E_tube[n], kdr_tube[n], NormalCrossEpo2[n],
//				Epo1[n], tPO1[n]);
//
//	}
//
//}
void calAllGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
		const Vector3 &Ein,const int &MaximumNumberOfReflections,
		 Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
		double *kdr_tube, int *hitnum_tube, Vector3 ** EpoA, Vector3 *Epo1,
		double *tPO1, double **kdotr, Vector3 **RefPoint){
	for (int n = 0; n < tnum; n++) {
		hitnum_tube[n] = GeometricalOptics(tube[n], bvh,  kin, Ein, MaximumNumberOfReflections,
				refpoint_tube[n], k_tube[n], E_tube[n], kdr_tube[n], EpoA[n],
				Epo1[n], tPO1[n], kdotr[n], RefPoint[n]);

	}

}






int main(int argc, char *argv[]) {
//	double Fmin = 14.5e9;
//	double Fmax = 15.5e9;
//	double Fstep = 0.01e9;
//	int NumberOfFrequencyPoints = (Fmax-Fmin)/Fstep+1;
//	cout<<NumberOfFrequencyPoints<<endl;
//	exit(0);

//-------------------------------input .3ds model---------------------------------
//	string input = "/home/user/cuda-workspace/lib3/Ball.3ds";
//	string input = "/home/user/cuda-workspace/lib3/plate1515.3ds";
	string input = "/home/user/cuda-workspace/lib3/dihedral 45d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral m.3ds";
//	string input = "/home/user/cuda-workspace/lib3/dihedral 22pt5d.3ds";
//	string input = "/home/user/cuda-workspace/lib3/trihedral30.3ds";
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

	cout<<"#Face ="<<num_face<<endl;
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



//	for (int ii =0 ; ii < NumberOfFrequencyPoints ; ii++){      //frequency loop
//		cout << (double)clock() / CLOCKS_PER_SEC << " S";

//		cout<<"n"<<ii<<endl;							 //frequency loop
//		double f;							//frequency loop
//		f = (ii*Fstep + Fmin);					 //frequency loop

		double f = 15e9;
		double lambda = c / f;
		cout<<"L="<<lambda<<endl;
//		exit(0);


		double Lx = 0.6;
		double Ly = 0.6;
		double rayPerWavelenth = 16;
		double cellSize = lambda/rayPerWavelenth;
		double Pmin = -45;
		double Pmax = 45;
		double Pstep = 0.5;
		int PointOfAngle = (Pmax-Pmin)/Pstep + 1;


		int xdnum = Lx / cellSize;
		int ydnum = Ly / cellSize;
		double k0 = 2 * PI * f * pow(c, -1);
		const int MaximumNumberOfReflections = 3; //Maximum hit-point depth
		double distance = 40;
//		distance > 2D^2/lambda
		string POL = "V";


		int rnum = xdnum * ydnum;
		int tnum = (xdnum - 1) * (ydnum - 1);

		Vector3* antr = new Vector3[rnum]; //ray  not used
		Vector3* tube = new Vector3[tnum]; //raytube
		Vector3* k_antr = new Vector3[rnum];//ray  not used
		Vector3* k_tube = new Vector3[tnum];
		Vector3* E_antr = new Vector3[rnum];//ray  not used
		Vector3* E_tube = new Vector3[tnum];
		Vector3* E_ii = new Vector3[tnum];
		Vector3* Epo1 = new Vector3[tnum];
		Vector3** EpoA = new Vector3*[tnum];
		  for(int i = 0; i < tnum; ++i){
			  EpoA[i] = new Vector3[MaximumNumberOfReflections];
		  }


		Vector3* kt_in = new Vector3[tnum];
		Vector3* Epo2 = new Vector3[tnum];
		Vector3* Epo3 = new Vector3[tnum];


		int *hitnum_antr = new int[rnum];//ray  not used
		int *hitnum_tube = new int[tnum];
		double* kdr_antr = new double[rnum]; // not used Point of emission to the point of reflection distance

		double* kdr_tube = new double[tnum];  //Point of emission to the point of reflection distance
		double** kdotr = new double*[tnum];  //Point of emission to the point of reflection distance
		  for(int i = 0; i < tnum; ++i){
			  kdotr[i] = new double[MaximumNumberOfReflections];
		  }
		double* kdotr2 = new double[tnum];  //Point of emission to the point of reflection distance
		double* kdotr3 = new double[tnum];



		double* kdr_antrr = new double[rnum]; //not used reflection point
		double* kdr_tuber = new double[tnum]; //reflection point to reflection plane distance
		double* total_path_tube = new double[tnum];
		double* kdrPO1= new double[tnum];

		double* IR_Angle = new double[tnum];


		CXV3* Ec_antr = new CXV3[rnum];//ray not used
		CXV3* Ec_tube = new CXV3[tnum];
		CXV3* EcPO1 = new CXV3[tnum];
//		CXV3* EcPO2 = new CXV3[tnum];
		CXV3* EcPO3 = new CXV3[tnum];



		double* RCSTheta = new double[PointOfAngle];
		double* RCSPhi =new double[PointOfAngle];
		Vector3* retpoint_antr = new Vector3[rnum]; //The Return point on incident plane
		Vector3* refpoint_antr = new Vector3[rnum]; //The source point on the object
		Vector3* retpoint_tube = new Vector3[tnum]; //The Return point on incident plane
		Vector3* refpoint_tube = new Vector3[tnum]; //The source point on the object
		Vector3** RefPoint = new Vector3*[tnum]; //The source point on the object
		  for(int i = 0; i < tnum; ++i){
			  RefPoint[i] = new Vector3[MaximumNumberOfReflections];
		  }
		Vector3* RefPoint2 = new Vector3[tnum];
		Vector3* RefPoint3 = new Vector3[tnum];


//		CXF* Atheta = new CXF[tnum];
//		CXF* Aphi = new CXF[tnum];
		CXF* Atheta_total = new CXF[MaximumNumberOfReflections+1];
		CXF* Aphi_total = new CXF[MaximumNumberOfReflections+1];

		int* OneAngleNumberOfHit = new int[PointOfAngle];
//		double* RCSABS = new double[PointOfAngle];
		CXF* RCSABS = new CXF[PointOfAngle];
		double* PHI = new  double[PointOfAngle];
		double* THETA = new double[PointOfAngle];
		double* ReE_Theta = new double[PointOfAngle];
		double* ImE_Theta = new double[PointOfAngle];
		double* ReE_Phi = new double[PointOfAngle];
		double* ImE_Phi = new double[PointOfAngle];

//

//		CXF* RCSABS = new CXF[phi2.size()];
//		double* PHI = new  double[phi2.size()];
//		double* THETA = new double[phi2.size()];
//		double* ReE_Theta = new double[phi2.size()];
//		double* ImE_Theta = new double[phi2.size()];
//		double* ReE_Phi = new double[phi2.size()];
//		double* ImE_Phi = new double[phi2.size()];
//		int* OneAngleNumberOfHit = new int[phi2.size()];
//		double* RCSTheta = new double[phi2.size()];
//		double* RCSPhi =new double[phi2.size()];




		bool *return_hit = new bool[tnum];

		memset(refpoint_antr, 0, rnum);
		memset(hitnum_tube, 0, tnum);
		memset(kt_in, 0, tnum);

		RayIndex *index;
		index = new RayIndex[tnum];




//		double PHI;
//		double THETA = -90 * PI / 180;
		Vector3 phi_h;
		Vector3 theta_h;
		Vector3 Ei;

//		for (int i = 0  ; i < phi2.size(); i++) {
		for (int i = 0  ; i < PointOfAngle; i++) {

			PHI[i] = (i*Pstep + Pmin) * PI / 180;
//			PHI[i] = phi2[i] * PI / 180;

			THETA[i] = -90 * PI / 180;

	//		PHI =  phi2[i]*PI/180;
			cout<<"PHI = "<<PHI[i]/PI*180<<endl;
//			cout<<"PHI = "<<PHI[i]<<endl;

			genrateAperture(THETA[i], PHI[i], distance, cellSize, rnum, tnum, Lx,
					Ly, xdnum, ydnum, antr, tube, kt_in[0]);



			initializeEF(POL, kt_in[0], theta_h, phi_h, Ei);
	//		cout<< "THETA_"<<"("<<theta_h.x<<","<<theta_h.y<<","<<theta_h.z<<")"<<endl;
	//		cout<< "PHI_"<<"("<<phi_h.x<<","<<phi_h.y<<","<<phi_h.z<<")"<<endl;
			calAllGeoOpt(tnum, tube ,bvh, kt_in[0], Ei, MaximumNumberOfReflections,
						refpoint_tube,  k_tube,  E_tube, kdr_tube, hitnum_tube, EpoA, Epo1, kdrPO1, kdotr,
						RefPoint);
//			exit(0);

	//
//			AllRayGeoOpt(tnum, tube ,bvh, kt_in[0], Ei, MaximumNumberOfReflections,
//					refpoint_tube,  k_tube,  E_tube, kdr_tube, hitnum_tube, Epo2, Epo1, kdrPO1);

			for(int iii=0 ;iii <tnum ; iii++){
				Epo1[iii] = EpoA[iii][0];
				Epo2[iii] = EpoA[iii][1];
				Epo3[iii] = EpoA[iii][2];
				kdrPO1[iii] = kdotr[iii][0];
				kdotr2[iii] = kdotr[iii][1];
				kdotr3[iii] = kdotr[iii][2];

				RefPoint2[iii] = RefPoint[iii][1];
				RefPoint3[iii] = RefPoint[iii][2];
			}

	// 	   convert to complex field
			AllEF(tnum, hitnum_tube, kt_in, k_tube, k_antr, index, xdnum,
					tube, RefPoint2, kdotr2,
					f, Epo2,
					kdr_tuber, retpoint_tube, total_path_tube,
					Ec_tube, IR_Angle);
			AllEF(tnum, hitnum_tube, kt_in, k_tube, k_antr, index, xdnum,
					tube, RefPoint3, kdr_tube,
					f, Epo3,
					kdr_tuber, retpoint_tube, total_path_tube,
					EcPO3, IR_Angle);



			PO1CEF(tnum, hitnum_tube, kt_in, k_tube, k_antr, index, xdnum,
					tube, refpoint_tube, kdrPO1,
					f, Epo1,
					kdr_tuber, retpoint_tube, total_path_tube,
					EcPO1, IR_Angle);



		calEquationTwentySeven(EcPO1, k0,
						cellSize, tnum, theta_h,
						phi_h, PHI[i], IR_Angle, retpoint_antr, Atheta_total[0], Aphi_total[0]);


		 calEquationTwentySeven(Ec_tube, k0,
					   cellSize, tnum, theta_h,
					   phi_h, PHI[i], IR_Angle, retpoint_antr, Atheta_total[1], Aphi_total[1]);

		 calEquationTwentySeven(EcPO3, k0,
						cellSize, tnum, theta_h,
						phi_h, PHI[i], IR_Angle, retpoint_antr, Atheta_total[2], Aphi_total[2]);

		 Atheta_total[3] = Atheta_total[0] + Atheta_total[1]+Atheta_total[2];
		 Aphi_total[3] = Aphi_total[0] + Aphi_total[1]+Aphi_total[2];

		 ReE_Theta[i] = Atheta_total[2].real();
		 ReE_Phi[i] = Aphi_total[2].real();
		 ImE_Phi[i] = Aphi_total[2].imag();
		 ImE_Theta[i] = Atheta_total[2].imag();




			RCSTheta[i] = 10* log10(PI * 4 * abs(Atheta_total[3]) * abs(Atheta_total[3]));
			RCSPhi[i]= 10 * log10(PI * 4 * abs(Aphi_total[3]) * abs(Aphi_total[3]));


			RCSABS[i] = 10 * log10(PI * 4 * abs(Aphi_total[3]) * abs(Aphi_total[3])+
					PI * 4 * abs(Atheta_total[3]) * abs(Atheta_total[3]));
//			RCSABS[i] =Aphi_total[1];
	//
			cout << "RCS(theta) ="
					<< RCSTheta[i]
					<< " dBsm" << endl;
			cout << "RCS(phi) ="<<
					RCSPhi[i]<< " dBsm"
					<< endl;
//			exit(0);




			getAngleData(PointOfAngle, RCSABS);

		}

	//	getAngleData(PointOfAngle, OneAngleNumberOfHit);
//	    ofstream outFile;
//	    int Number_of_files=PointOfFrequency;
//	    string CSTformat;
//
//        CSTformat="rcs" + IntToStr(ii) +".dat";
//        cout<< CSTformat << "  \n";
//		outFile.open(CSTformat.c_str());
//        for (int i = 0; i < PointOfAngle; i++) {
//
//			outFile << RCSPhi[i]<<endl;
//        }
//			outFile.close();
//		getCSTformatData(THETA, PHI, ReE_Theta, ImE_Theta,
//				ReE_Phi, ImE_Phi, PointOfAngle, NumberOfFrequencyPoints, ii);
//		getCSTformatData(THETA, PHI, ReE_Theta, ImE_Theta,
//						ReE_Phi, ImE_Phi, phi2.size(), NumberOfFrequencyPoints, ii);


		getMatlabData(POL, PointOfAngle, RCSTheta, RCSPhi);

		delete[] antr;
		delete[] tube;
		delete[] k_antr;
		delete[] k_tube;
		delete[] E_antr;
		delete[] E_tube;
		delete[] hitnum_antr;
		delete[] hitnum_tube;
		delete[] kdr_antr; //distance Point of emission to the point of reflection
		delete[] kdr_tube; //Point of emission to the point of reflection
		delete[] kdr_antrr; //reflection point
		delete[] kdr_tuber; //reflection point

		delete[] Ec_antr;
		delete[] Ec_tube;
		delete[] Epo1;
		delete[] EcPO1;


		for(int i = 0; i<tnum ; ++i){

		    delete [] EpoA[i];
		  }
		delete[] EpoA;
		delete[] Epo3;


		for(int i = 0; i<tnum ; ++i){

		    delete [] kdotr[i];
		  }
		delete[] kdotr;
		delete[] kdotr2;
		delete[] kdotr3;

		delete[] retpoint_antr; //The Return point on incident plane
		delete[] refpoint_antr; //The source point on the object
		delete[] retpoint_tube; //The Return point on incident plane
		delete[] refpoint_tube; //The source point on the object

		for(int i = 0; i<tnum ; ++i){

		    delete [] RefPoint[i];
		  }
		delete [] RefPoint;
		delete [] RefPoint2;
		delete [] RefPoint3;
//		delete[] Atheta;
//		delete[] Aphi;
		delete[] total_path_tube;
		delete[] return_hit;
		delete[] index;
		delete[] kt_in;
		delete[] Atheta_total;
		delete[] Aphi_total;
		delete[] RCSTheta;
		delete[] RCSPhi;
		delete[] IR_Angle;
		delete[] OneAngleNumberOfHit;
		delete[] E_ii;
		delete[] Epo2;
		delete[] EcPO3;
		delete[] RCSABS;
		delete[] kdrPO1;
		delete[] PHI;
		delete[] THETA;
		delete[] ReE_Theta;
		delete[] ImE_Theta;
		delete[] ReE_Phi;
        delete[] ImE_Phi;

	cout << (double)clock() / CLOCKS_PER_SEC << " S";

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
void ComplexField(const Vector3 &R, const double &kr, const double frequency, CXV3 &C) {
	double multiply2 = 2.0;
	CXF expvalue = exp(j * (kr) / c * frequency * multiply2 * PI);
	C = CXV3(CXF(R.x, 0) * expvalue, CXF(R.y, 0) * expvalue,
			CXF(R.z, 0) * expvalue);
//	C = CXV3(CXF(R.x, 0) , CXF(R.y, 0) , CXF(R.z, 0) );

}

double Vector_Angle_d(const Vector3 &v1, const Vector3 &v2) {
	double l1 = length(v1);
	double l2 = length(v2);
	double angle = 0;
//	cout << v2.x << "," << v2.y << "," << v2.z << endl;

//	cout<<acos(-1) <<endl;

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
	delete[] R2;
}




//void AllRayGeoOpt(const int &tnum, const Vector3 *tube, const BVH &bvh, const Vector3 &kin,
//		const Vector3 &Ein,const int &MaximumNumberOfReflections,
//		 Vector3 *refpoint_tube, Vector3 *k_tube, Vector3 *E_tube,
//		double *kdr_tube, int *hitnum_tube, Vector3 *NormalCrossEpo2, Vector3 *Epo1,
//		double *tPO1){
//	for (int n = 0; n < tnum; n++) {
//		hitnum_tube[n] = GeometricalOptics(tube[n], bvh,  kin, Ein, MaximumNumberOfReflections,
//				refpoint_tube[n], k_tube[n], E_tube[n], kdr_tube[n], NormalCrossEpo2[n],
//				Epo1[n], tPO1[n]);
//
//	}
//
//}


void calEquationTwentySeven(const CXV3 *E, const double &k0,
		const double &cellSize, const int &tnum, const Vector3 &theta_h,
		const Vector3 &phi_h, const double &PHI, const double *IR_Angle, Vector3 *retpoint_antr,
		CXF &Atheta, CXF &Aphi){
	CXF* SingleRayAtheta = new CXF[tnum];
	CXF* SingleRayAphi = new CXF[tnum];
	CXF expvalue = CXF(1,0);
	CXF Ex, Ey;
	double Ii = 1.0;
	//	float Sx = 0;
	//	float Sy = 0;
	//	float Xi = retpoint * phi_h;
	//	float Yi = retpoint * theta_h;
	//	cout<<Xi<<","<<Yi<<endl;
	//	cout << "E=" << E.x << "," << E.y << "," << E.z << endl;
	//	CXF expvalue = exp(j * k0 * (Sx * Xi * +Sy * Yi));

	Atheta = CXF(0, 0);
	Aphi = CXF(0, 0);
	for(int i = 0 ; i < tnum ; i++){

			Ex = E[i] * phi_h;
//			if ( Ex!=CXF(0,0)){
//			cout<<"Ex ="<<Ex<<endl;
//			}
			Ey = E[i] * theta_h;

//			if ( Ey!=CXF(0,0)){
//		 	cout<<"Ey ="<<Ey<<endl;
//			}

		//	float Phi_p; //Phi for paper

			if (Ex == CXF(0, 0) && Ey == CXF(0, 0)) {

				SingleRayAtheta[i] = CXF(0, 0);
				SingleRayAphi[i] = CXF(0, 0);
			} else {

		//		SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos((180-IR_Angle[i])*PI/180);
		//		SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos(PHI)/cos(PHI);
		//		SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)/cos(PHI);

		//	   SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii)*(1+tan(PHI)*tan(2*PHI));
			   SingleRayAtheta[i] = j*k0/(2*PI)*( Ey ) * expvalue * cellSize * cellSize* (Ii);

		//			cout<<"Atheta= "<<Atheta<<endl;
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos((180-IR_Angle[i])*PI/180);
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(PHI)/cos(PHI);
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(PHI);

		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)/cos(IR_Angle[i]*PI/180);
		//		SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii)*(1+tan(PHI)*tan(2*PHI));
				SingleRayAphi[i] = j*k0/(2*PI)*(Ex ) * expvalue * cellSize * cellSize * (Ii);



		//			cout<<"Aphi= "<<Aphi<<endl;

			}



			Atheta += SingleRayAtheta[i];
//			cout<<"single = "<<SingleRayAtheta[i]<<endl;
//			cout<<Atheta<<endl;
			Aphi += SingleRayAphi[i];
//			cout<<Aphi<<endl;

	}
	delete[] SingleRayAtheta;
	delete[] SingleRayAphi;


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

void genrateAperture(const double&THETA, const double &PHI, const double &distance,
		const double cellSize, const int rnum,const int tnum,
		const double &Lx, const double &Ly, const int xdnum, const int ydnum, Vector3 *antr, Vector3 *antt,
		Vector3  &ant_N){
	Vector3 antc, cros_1, cros_2, ant_Nd;

	for (int i = 0; i < rnum; i++) {
		antr[i].x = 0 + (i % xdnum) * cellSize;
		antr[i].y = 0 + i / xdnum * cellSize;
		antr[i].z = 0;
	}

	for (int i = 0; i < tnum; i++) {
		antt[i].x = (0 + (i % (xdnum - 1)) * cellSize + cellSize / 2);
		antt[i].y = (0 + i / (xdnum - 1) * cellSize + cellSize / 2);
		antt[i].z = 0;
	}
//
//	for (int i = 0; i < rnum; i++) {
//		  cout<<antr[i].x<<","<<antr[i].y<<","<<antr[i].z<<"antr"<<endl;
//	}

	antc = (antr[0] + antr[rnum - 1]) / 2.0;


	for (int i = 0; i < rnum; i++) {
		antr[i] = antr[i] - antc;
	}


	for (int i = 0; i < tnum; i++) {
		antt[i] = antt[i] - antc;
	}

	rotate_THETA_then_PHI(antr, rnum, THETA, PHI);
	rotate_THETA_then_PHI(antt, tnum, THETA, PHI);


	cros_1 = antr[1] - antr[0];
	cros_2 = antr[rnum] - antr[0];

	ant_N = cros_2 ^ cros_1;
	ant_N = normalize(ant_N);
	ant_Nd = distance * ant_N;

	for (int i = 0; i < rnum; i++) {
		antr[i] = antr[i] - ant_Nd;
	}

	for (int i = 0; i < tnum; i++) {

		antt[i] = antt[i] - ant_Nd;
	}


//	antc = antc - ant_Nd;

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

void AllEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
		const Vector3 *tube, const Vector3 *refpoint_tube, const double *kdr_tube,
		const double f, const Vector3 *E_tube,
		double *kdr_tuber, Vector3 *retpoint_tube, double *total_path_tube,
		CXV3 *Ec_tube, double *IR_Angle){

//	int a, b;

	for (int n = 0; n < tnum; n++) {
//		a = n / (xdnum - 1);
//		b = n % (xdnum - 1);
//
//		index->ray_center = n;
//		(*index).ray_one = xdnum * a + b;
//		(*index).ray_two = (*index).ray_one + 1;
//		(*index).ray_three = (*index).ray_two + xdnum;
//		(*index).ray_four = (*index).ray_one + xdnum;


		if(hitnum_tube[n]!=0){
			IR_Angle[n] = Vector_Angle_d(kt_in[0], k_tube[n]);
//				cout<<IR_Angle[n]<<endl;

		}else{
			IR_Angle[n] = -1;
		}
		if (IR_Angle[n] > 0 ){
			intersectPlane((-1 * kt_in[0]), tube[0],
					refpoint_tube[n], (-1 * kt_in[0]), kdr_tuber[n]);
//			intersectPlane((-1 * kt_in[0]), tube[0],
//					refpoint_tube[n], k_tube[n], kdr_tuber[n]);

			retpoint_tube[n] = refpoint_tube[n] + kdr_tuber[n] * k_tube[n];
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
void PO1CEF(const int &tnum, const int *hitnum_tube,const Vector3 *kt_in, const Vector3 *k_tube, Vector3 *k_antr,
		RayIndex *index,const int xdnum,
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


//void Raytube_Numbering(const int &tnum, const int &xdnum, const Vector3 *antr, const Vector3 *antt, RayIndex *index) {
//
////	RayIndex *index;
////	index = new RayIndex[tnum];
//
//	int m, n;
//	for (int i = 0; i < tnum; i++) {
//		m = i / (xdnum - 1);
//		n = i % (xdnum - 1);
//
//		index->ray_center = i;
//		(*index).ray_one = xdnum * m + n;
//		(*index).ray_two = (*index).ray_one + 1;
//		(*index).ray_three = (*index).ray_two + xdnum;
//		(*index).ray_four = (*index).ray_one + xdnum;
//
////
////    cout << index->ray_center<<","\
////         << (*index).ray_one << ","\
////         << (*index).ray_two << ","\
////         << (*index).ray_three << ","\
////         << (*index).ray_four << endl;
////
////
////
////  cout << antt[index->ray_center].x << ", " << antt[index->ray_center].y << "," << antt[index->ray_center].z  << "tube" << endl;
////  cout << antr[(*index).ray_one].x << ", " << antr[(*index).ray_one].y << "," << antr[(*index).ray_one].z  << "antr_one" << endl;
////  cout << antr[(*index).ray_two].x << ", " << antr[(*index).ray_two].y << "," << antr[(*index).ray_two].z  << "antr_two" << endl;
////  cout << antr[(*index).ray_three].x << ", " << antr[(*index).ray_three].y << "," << antr[(*index).ray_three].z  << "antr_three" << endl;
////  cout << antr[(*index).ray_four].x << ", " << antr[(*index).ray_four].y << "," << antr[(*index).ray_four].z  << "antr_four" << endl;
//
//	}
//}




