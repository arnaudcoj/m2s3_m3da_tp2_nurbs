#include "Nurbs.h"
#include <cmath>
#include <iostream>
#include <algorithm>

#include "Vector3.h"
#include "Matrix4.h"

using namespace std;
using namespace p3d;

Nurbs::~Nurbs() {
}

Nurbs::Nurbs() {
  _control.clear();
  _nbControl[D_U]=0;
  _nbControl[D_V]=0;
  _knot[D_U].clear();
  _knot[D_V].clear();
  _degree[D_U]=0;
  _degree[D_V]=0;

}

double Nurbs::startInterval(EDirection direction) {
  return _knot[direction][_degree[direction]];
}

double Nurbs::endInterval(EDirection direction) {
  return _knot[direction][_nbControl[direction]]-0.00001;
}


bool Nurbs::inInterval(EDirection direction, double u) {
  return (u>=startInterval(direction) && u<=endInterval(direction));
}


void Nurbs::knotUniform(EDirection direction,int nb) {
  _knot[direction].resize(nb);
  /* TODO : set uniform knots
   * _knot[direction][i] is the knot u_i for the given direction
   *
   *
   */
  double step = 1. / double(nb -1); // -1 car permet d'avoir la valeur 1 dans le tableau
  double cpt = 0.;
  for(int i = 0; i < nb; i++) {
    _knot[direction][i] = cpt;
    cpt += step;
  }

}


/** Eval the basis function Nkp(t) for the knot vector knot **/
double Nurbs::evalNkp(int k,int p,double u,std::vector<double> &knot) {
  double result=0.0;

  if(p == 0) {

    if(u >= knot[k] && u < knot[k + 1]) {
      result = 1.;
    } else {
      result = 0;
    }
  } else {
    double uk = knot[k];
    double ukp = knot[k + p];
    double upk1 = knot[p + k + 1];
    double uk1 = knot[k + 1];

    if((ukp - uk) != 0.)
      result += (u - uk) / (ukp - uk) * evalNkp(k, p -1, u, knot);
    if((upk1 - uk1) != 0.)
      result += (upk1 - u) / (upk1 - uk1) * evalNkp(k + 1, p - 1, u, knot);

  }
  return result;
}


double Nurbs::evalNkp(EDirection direction,int k,int p,double t) {
  return evalNkp(k,p,t,_knot[direction]);
}


void Nurbs::clearControl() {
  _nbControl[D_U]=0;
  _nbControl[D_V]=0;
  _control.clear();
}

void Nurbs::initControlGrid() {
  _nbControl[D_U]=5;
  _nbControl[D_V]=4;
  _control.clear();
  double u=-1;
  double v=-1;
  double stepU=2.0/(_nbControl[D_U]-1);
  double stepV=2.0/(_nbControl[D_V]-1);
  for(int i=0;i<_nbControl[D_V];++i) {
    u=-1;
    for(int j=0;j<_nbControl[D_U];++j) {
      _control.push_back(Vector4(u,v,double(rand())/RAND_MAX-0.5,1));
      u+=stepU;
    }
    v+=stepV;
  }
  knotRemap(D_U);
  knotRemap(D_V);
}


void Nurbs::addControlU(const Vector4 &p) {
  _control.push_back(p);
  _nbControl[D_U]++;
  knotRemap(D_U);
}


Vector3 Nurbs::pointCurve(double u) {
  Vector4 result(0,0,0,0);

  int n = nbControl(D_U);
  int p = degree(D_U);

  for(int k = 0; k < n; k++) {
    Vector4 pk = control(k);
    result += evalNkp(D_U, k, p, u) * pk;
  }

  return Vector3(result.x(),result.y(),result.z()) / result.w();
}


Vector3 Nurbs::pointSurface(double u,double v) {
  Vector4 result(0,0,0,0);

  return result.project(); // divide by w
}



void Nurbs::knotRemap(EDirection direction) {
  while (!checkNbKnot(direction)) {
    int nb=nbKnot(direction);
    _knot[direction].push_back(_knot[direction][nbKnot(direction)-1]);
    for(unsigned int i=nb-1;i>0;--i) {
      _knot[direction][i]=_knot[direction][i+1]-(_knot[direction][i]-_knot[direction][i-1])*(nb-1)/nb;
    }
  }
}


bool Nurbs::checkNbKnot(EDirection direction) {
  return (nbKnot(direction)>=nbControl(direction)+degree(direction)+1);
}


void Nurbs::knotOpenUniform(EDirection direction) {
  int nb = nbControl(direction);
  int deg = degree(direction);

  double step = 1. / (nb - deg);

  _knot[direction].resize(nb + deg + 1);

  int k = 0;

  for(int i = 0; i < deg + 1; i++, k++) {
    _knot[direction][k] = 0.;
  }

  for(int i = 1; i < nb - deg ; i++, k++) {
    _knot[direction][k] = i * step;
  }

  for(int i = 0; i < deg + 1; i++, k++) {
    _knot[direction][k] = 1.;
  }

  /*cout << nb << " " << deg << " " << k - 1 << " " << step << endl;
  for(k = 0; k < nb + deg + 1; k++) {
    cout << _knot[direction][k] << " ";
  }
  cout << endl;
  */
}


void Nurbs::knotBezier(EDirection direction) {
  int nb = nbControl(direction);
  degree(direction, nb - 1);

  int deg = degree(direction);

  _knot[direction].resize(nb + deg + 1);

  int k = 0;

  for(int i = 0; i < deg +1; i++, k++) {
    _knot[direction][k] = 0.;
  }
  for(int i = 0; i < deg +1; i++, k++) {
    _knot[direction][k] = 1.;
  }

  /*cout << nb << " " << deg << " " << k << " " << endl;
  for(k = 0; k < nb + deg + 1; k++) {
    cout << _knot[direction][k] << " ";
  }
  cout << endl;*/
}

void Nurbs::setCircle() {
  /* Have to set : _control, _degree[D_U], _knot[D_U], _nbControl[D_U]
   *
   */
  _control.clear();

}


void Nurbs::setRevolution(int nbV) {
  if (nbV<2) return;
  _nbControl[D_V]=nbV;
  _degree[D_V]=_degree[D_U];
  _control.resize(_nbControl[D_U]*_nbControl[D_V]);
  knotRemap(D_V);
  knotOpenUniform(D_V);

  double stepTheta=360.0/(nbV-1);
  double theta=stepTheta;
  Matrix4 rotate;
  for(int slice=nbControl(D_U);slice<nbControl(D_U)*nbControl(D_V);slice+=nbControl(D_U)) {
    rotate.setRotation(theta,0,1,0);
    for(int istack=0;istack<nbControl(D_U);++istack) {
      _control[slice+istack]=rotate*_control[istack];
    }
    theta+=stepTheta;
  }
}



