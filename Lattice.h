/* 
 * File:   Lattice.h
 * Author: qfeuille
 *
 * Created on 12 September 2012, 22:55
 */

#ifndef LATTICE_H
#define	LATTICE_H

#include <Eigen/Dense>

namespace NumMethod{
template<typename NetworkT, int Ndim>
class RectangularLattice {
    typedef Eigen::Matrix<int, Ndim, 1> Vec;
    
    public:
        NetworkT network;
        Vec dimensions;
        int getn (Vec position);
        Vec getPosition (int n);
        int getCoordinate (int n, int coord);
        int size (){return dimensions.prod();};
        int wrap (int x, int width);
        void operator () (int in, std::vector<int>& out);
        RectangularLattice (const Vec& idimensions);
      
};

template<typename NetworkT, int Ndim>
inline int RectangularLattice<NetworkT,Ndim>::getn (Vec position){
    int n =position [Ndim-1];
    for (int i=Ndim-1; i>0;i--){
        n = (n*dimensions [i-1])+position [i-1];
    }
    return n;
}

template<typename NetworkT, int Ndim>
inline Eigen::Matrix<int,Ndim,1> RectangularLattice<NetworkT,Ndim>::getPosition (int n){
    Vec position;
    double divfactor =1.0;
    int modfactor=1;
    for (int i=0; i<Ndim-1;i++){
        modfactor *= dimensions [i];
        position [i] = (n % modfactor)/divfactor; 
        n -= position [i];
        divfactor = (double) modfactor;
    }
    position [Ndim-1]= n/divfactor;
    return position;
    
}

//Do not use until fixed!
template<typename NetworkT, int Ndim>
inline int RectangularLattice<NetworkT,Ndim>:: RectangularLattice::getCoordinate (int n, int coord){
    int position;
    int divfactor =1;
    for (int i=0; i<coord;i++){
        position = (n % dimensions [i])/divfactor;
        n -= position;
        divfactor *= dimensions [i];
    }
    return position;
    
}

template<typename NetworkT, int Ndim>
void RectangularLattice<NetworkT,Ndim>::operator ()(int in, std::vector<int>& out){
    Vec position = getPosition (in);
        for (int i = 0; i < Ndim; i++) {
            int current = position [i];
            position[i] = wrap(current + 1, dimensions[i]);
            out.push_back( getn(position));
            position[i] = wrap(current - 1, dimensions[i]);
            out.push_back( getn(position));
            position [i] = current;
          
        } 
}

template<typename NetworkT, int Ndim>
inline int RectangularLattice<NetworkT,Ndim>::wrap (int x, int width){
    return (x%width + width)%width;
}

template<typename NetworkT, int Ndim>
RectangularLattice<NetworkT,Ndim>::RectangularLattice(const Vec& idimensions):
        network (){
    dimensions=idimensions;
    network.addNodes (idimensions.prod());
    network.setNodeConnections (*this);
}




}
#endif	/* LATTICE_H */

