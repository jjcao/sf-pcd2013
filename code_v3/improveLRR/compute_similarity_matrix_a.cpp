//mex -I"D:\PCL1.6.0\3rdParty\Eigen\include" -g compute_similarity_matrix_a.cpp

#include "mex.h" 
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <vector>
#include <algorithm> 
#include <iostream>
using namespace Eigen;
using namespace std;
#pragma warning( disable : 4267 24) 

// input argument 0: Xid
int nSample; double* Xid;
// input argument 1: neighbors(numberall*nSample)
double *neighbors;int numberall;
// input argument 2: Params
int npts;// number of points, all points
double *origin,*sigms;
double *Allpoints;//3*npts
// input argument 3: TP. TP.k_knn_normals is numberall
int numberi; double guide_rate,tolerate_angle,feature_threshold;
// output: similarity matrix
double *simiMat;

struct myclass {
	bool operator() (int i,int j) { return (i<j);}
} myobject;
void compute_svd()
{
	MatrixXd nhPts(3,numberall);// neighbors of a sample	
	double* nhData = nhPts.data();

	int i,j, t1,t2,t3;
	double tmp;
	int actualJ;
	//double* base = new double(3*nSample);
	for(i = 0; i<nSample; ++i)
	{
		t1 = i*numberall;
		double* nhi = &neighbors[t1];//neighbors(:,i)
		
		actualJ = 0;
		for(j = 0; j < numberall; ++j,++actualJ)
		{
			t2 = actualJ*3; t3 = (nhi[j]-1)*3;
			nhData[t2] = Allpoints[t3]-origin[0];
			nhData[t2+1]=Allpoints[t3+1]-origin[1];
			nhData[t2+2]=Allpoints[t3+2]-origin[2];
			tmp = sqrt(nhData[t2]*nhData[t2]+
				nhData[t2+1]*nhData[t2+1]+nhData[t2+2]*nhData[t2+2]);
			if(tmp==0)
				--actualJ;
			else
			{
				nhData[t2] = nhData[t2]/tmp;
				nhData[t2+1] = nhData[t2+1]/tmp;
				nhData[t2+2] = nhData[t2+2]/tmp;
			}
		}
		MatrixXd m = nhPts.block(0,0,3,actualJ);//neighbor points after filtering
		JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
		// in matlab, the default is from large to small, the same here
		//svd.singularValues();
		//cout << " left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
		//cout << " right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
		Vector3d v = svd.matrixU().col(2);//base(:,ri) = U1(:,size(tX,1));
		//mexPrintf("%f, %f, %f\n", v[0],v[1],v[2]);

	}// end of for(i = 0; i<nSample; ++i)
}
//estimation of normal of one point with randomn-neighbours
void computePointNormal(vector<int>indices,Vector3d &pointNormal)
{
	MatrixXd nhPts(3,indices.size());// neighbors of a sample	
	double* nhData = nhPts.data();
	int i,t1,t2;
	int actualJ=0;
	for(i = 0; i<indices.size(); ++i,++actualJ)
	{
		t1 = (indices[i]-1)*3;//标号是从1开始的
		t2 = actualJ*3; 
		nhData[t2] = Allpoints[t1]-origin[0];
		nhData[t2+1]=Allpoints[t1+1]-origin[1];
		nhData[t2+2]=Allpoints[t1+2]-origin[2];
	}
		MatrixXd m1 = nhPts.block(0,0,3,actualJ);//neighbor points after filtering
		MatrixXd m2 = m1 * m1.transpose();
		SelfAdjointEigenSolver<Matrix3d> eigensolver(m2);
		pointNormal = eigensolver.eigenvectors().col(0);//base(:,ri) = U1(:,size(tX,1));
}
//estimation of normal of sample point with neighbours
void compute_eig(vector<Vector3d>& normals)
{
	MatrixXd nhPts(3,numberall);// neighbors of a sample	
	double* nhData = nhPts.data();
	normals.resize(nSample);
	int i,j, t1,t2,t3;
	double tmp;
	int actualJ;
	//double* base = new double(3*nSample);
	for(i = 0; i<nSample; ++i)
	{
		t1 = i*numberall;
		double* nhi = &neighbors[t1];//neighbors(:,i)
		
		actualJ = 0;
		for(j = 0; j < numberall; ++j,++actualJ)
		{
			t2 = actualJ*3; t3 = (nhi[j]-1)*3;
			nhData[t2] = Allpoints[t3]-origin[0];
			nhData[t2+1]=Allpoints[t3+1]-origin[1];
			nhData[t2+2]=Allpoints[t3+2]-origin[2];
			/*tmp = sqrt(nhData[t2]*nhData[t2]+
			nhData[t2+1]*nhData[t2+1]+nhData[t2+2]*nhData[t2+2]);
			if(tmp==0)
			--actualJ;
			else
			{
			nhData[t2] = nhData[t2]/tmp;
			nhData[t2+1] = nhData[t2+1]/tmp;
			nhData[t2+2] = nhData[t2+2]/tmp;
			}*/
		}
		MatrixXd m1 = nhPts.block(0,0,3,actualJ);//neighbor points after filtering
		MatrixXd m2 = m1 * m1.transpose();
		SelfAdjointEigenSolver<Matrix3d> eigensolver(m2);
		// in matlab, the default is from small to large, the same here
//cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << endl;
//cout << "Here's a matrix whose columns are eigenvectors of A \n"
//<< "corresponding to these eigenvalues:\n"
//<< eigensolver.eigenvectors() << endl;

		Vector3d v0 = eigensolver.eigenvectors().col(0);//base(:,ri) = U1(:,size(tX,1));
		normals[i]=v0;
		/*float e0=eigensolver.eigenvalues()[0];
		Vector3d v1 = eigensolver.eigenvectors().col(1);
		float e1=eigensolver.eigenvalues()[1];
		Vector3d v2 = eigensolver.eigenvectors().col(2);
		float e2=eigensolver.eigenvalues()[2];*/
		//mexPrintf("1th eigenvector=%f, %f, %f,1th eigenvalue=%f\n",v0[0],v0[1],v0[2],e0);
		//mexPrintf("2th eigenvector=%f, %f, %f,2th eigenvalue=%f\n",v0[0],v0[1],v0[2],e1);
		//mexPrintf("3th eigenvector=%f, %f, %f,3th eigenvalue=%f\n",v0[0],v0[1],v0[2],e2);
	}// end of for(i = 0; i<nSample; ++i)
}
void  ransacNumber(vector<int>&rnumber,double *neigh_one)
{
	rnumber.resize(numberi);
	//mexPrintf("numberi=%d\n",numberi);
	vector<int>flag(numberall,0);
	vector<double> p(numberall);double sum=0.0;
	int cont=0;//record choosed point
	for (int i=0;i<numberall;i++)
	{
		int a=neigh_one[i];
		//mexPrintf("a=%d\n",a);
		double a_sigma=exp(-sigms[a-1]);sum+=a_sigma;
		p[i]=sum;
		//mexPrintf("p(%d)=%f\n",i,p[i]);
	}
	while (cont<numberi)
	{
		double s=RAND_MAX;
		double v=sum*(rand()/s);
		//double v=sum*((rand()%1000)/1000.0);
		//mexPrintf("v=%d\n",v);
		for (int i=0;i<numberall;i++)
		{
			if (p[i]>=v)
			{
				if (flag[i]==1)
				{
					break;
				}
				flag[i]=1;
				//int bk=neigh_one[i];
				rnumber[cont]=neigh_one[i];
				//mexPrintf("rand_one: %d\n",rnumber[cont]);
				cont++;
				break;
			} 
		}
	}
}
void  compute_similarity_matrix()
{
	int ti,tj,tm;double m1,m2;
	Vector3d normal1,normal2;
	//MatrixXd simMatrix(nSample,nSample);
	//vector<double>sim_value(nSample);
	vector<double>sim_value2(nSample);
	vector<Vector3d>normals;
	compute_eig(normals);
	for (int i=0;i<nSample;i++)
	{
		ti=i*numberall;
		tm=Xid[i]-1;
		m1=sigms[tm];
		double*nthi=&neighbors[ti];
		//simMatrix(i,i)=1.0;
		simiMat[i*nSample+i]=1.0;
		for (int j=i+1;j<nSample;j++)
		{
			tj=j*numberall;
			tm=Xid[j]-1;
			m2=sigms[tm];
			double* nthj=&neighbors[tj];
			//计算随机标号
			if (m1>feature_threshold)
			{
				vector<int>rnumber_i(numberi);
				ransacNumber(rnumber_i,nthi);
				computePointNormal(rnumber_i,normal1);
			}
			else
			{
				normal1=normals[i];
			}
			if (m2>feature_threshold)
			{
				vector<int>rnumber_j(numberi);
				ransacNumber(rnumber_j,nthj);
				computePointNormal(rnumber_j,normal2);
			}
			else
			{
				normal2=normals[j];
			}
			//simMatrix(i,j)=abs((double)normal1.dot(normal2));
			//simMatrix(j,i)=simMatrix(i,j);
			simiMat[i*nSample+j]=abs((double)normal1.dot(normal2));
			simiMat[j*nSample+i]=simiMat[i*nSample+j];
			//test
			/*for (int k1=0;k1<numberall;k1++)
			{
				mexPrintf("all:%f\n",nthi[k1]);
				int ll=nthi[k1];
				mexPrintf("sigma:%f\n",sigms[ll-1]);
			}
			for (int kk=0;kk<rnumber_i.size();kk++)
			{
				mexPrintf("rand:%d\n",rnumber_i[kk]);
			}*/
		}
		/*for (int e1=0;e1<nSample;e1++)
		{
			sim_value[e1]=simMatrix(i,e1);
		}*/
		for (int e1=0;e1<nSample;e1++)
		{
			sim_value2[e1]=simiMat[i*nSample+e1];
		}
		//sort(sim_value.begin(),sim_value.end(),myobject);
		sort(sim_value2.begin(),sim_value2.end(),myobject);
		//double throled=max(sim_value[(int)nSample*guide_rate],tolerate_angle);
		double throled2=max(sim_value2[(int)nSample*guide_rate],tolerate_angle);
	/*	for (int e2=0;e2<nSample;e2++)
		{
			if (simMatrix(i,e2)<throled)
			{
				simMatrix(i,e2)=1.0;
			} 
			else
			{
				simMatrix(i,e2)=0.0;
			}
		}*/
		for (int e2=0;e2<nSample;e2++)
		{
			if (simiMat[i*nSample+e2]<throled2)
			{
				simiMat[i*nSample+e2]=1.0;
			} 
			else
			{
				simiMat[i*nSample+e2]=0.0;
			}
		}
	}
	//test show simMatrix
	/*for (int i=0;i<nSample;i++)
	{
		for(int j=0;j<nSample;j++)
		{
			mexPrintf("%f ",simMatrix(j,i));
		}
		mexPrintf("\n");
	}*/
}
void mexFunction(int nlhs, mxArray *plhs[],	int nrhs, const mxArray *prhs[]) 
{
	///////////// Error Check
	if (4 != nrhs) 
		mexErrMsgTxt("Number of input should be 4");

	if (1 != nlhs) 
		mexErrMsgTxt("Number of output should be 1");


	///////////// input & output arguments	
	// input 0: Xid: 1*nSample
	nSample = mxGetN(prhs[0]);
	Xid = mxGetPr(prhs[0]); // retrieve input
	
	// output 0
	plhs[0] = mxCreateDoubleMatrix(nSample,nSample,mxREAL);// output		
	simiMat = mxGetPr(plhs[0]);

	// input 1: neighbours: numberall*nSample
	numberall = mxGetM(prhs[1]);
	neighbors = mxGetPr(prhs[1]); 
	//mexPrintf("%f, %f, %f\n", neighbors[0],neighbors[1],neighbors[2]);

	// input 2: Params	
	mxArray* tmp;
	const mxArray *Params = prhs[2];
	if ( mxSTRUCT_CLASS != mxGetClassID(Params))
		mexErrMsgTxt("3rd arguments is not a structure!");
	else
	{
		// Params.origin: 3*1
		tmp = mxGetField(Params,0,"origin");
		origin = mxGetPr(tmp);
		//mexPrintf("%f, %f, %f\n", origin[0],origin[1],origin[2]);

		// Params.sigms: npts*1
		tmp = mxGetField(Params,0,"sigms");
		npts = mxGetM(tmp);
		sigms = mxGetPr(tmp);

		// Params.Allpoints: 3*npts
		tmp = mxGetField(Params,0,"Allpoints");
		Allpoints = mxGetPr(tmp);		
	}

	// input 3: TP
	const mxArray *TP = prhs[3];
	if ( mxSTRUCT_CLASS != mxGetClassID(TP))
		mexErrMsgTxt("3rd arguments is not a structure!");
	else
	{
		// numberi = TP.numberi;
		tmp = mxGetField(TP,0,"numberi");
		numberi = mxGetScalar(tmp);
		//mexPrintf("%d\n", numberi);

		// guide_rate = TP.guide_rate;
		tmp = mxGetField(TP,0,"guide_rate");
		guide_rate = mxGetScalar(tmp);

		// tolerate_angle = TP.tolerate_angle;
		tmp = mxGetField(TP,0,"tolerate_angle");
		tolerate_angle = mxGetScalar(tmp);

		//TP.feature_threshold
		tmp = mxGetField(TP,0,"feature_threshold");
		feature_threshold = mxGetScalar(tmp);
	}

	//compute_svd();
	/*vector<Vector3d>normals;Vector3d a;
	compute_eig(normals);
	for (int i=0;i<normals.size();i++)
	{
	mexPrintf("%d th eigenvector=%f, %f, %f\n",i,normals[i](0),normals[i](1),normals[i](2));
	vector<int>indices(numberall);
	double *k=&neighbors[i*numberall];
	for (int j=0;j<numberall;j++)
	{
	indices[j]=k[j];
	}
	computePointNormal(indices,a);
	mexPrintf("%d th eigenvector test_computePointNormal=%f, %f, %f\n",i,a(0),a(1),a(2));
	}*/
	 compute_similarity_matrix();
}
