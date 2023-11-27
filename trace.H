#ifndef TRACE_H
#define TRACE_H

#include "slVector.H"
#include <vector>




class Ray {  //光线向量
public:
    SlVector3 e;
    SlVector3 d;
    int depth;
    Ray(const SlVector3 &_e, const SlVector3 &_d, int _depth = 0) : e(_e), d(_d), depth(_depth) {};
};

class Fill {  //填充
public: 
    SlVector3 color;
    double kd, ks, shine, t, ior;
};

class HitRecord {  //光照到的东西
public:
    double t, alpha, beta, gamma;
    SlVector3 p, n, e, v;
    Fill f;
    int raydepth;
    SlVector3  getRefractionD() const;
};

class Light {
public:
    SlVector3 p, c;
};

struct Node;
class Box {
public:
    double bx, sx, by, sy, bz, sz;
    std::vector<int> intersectList;
    Box() {}
    Box(double _bx, double _sx, double _by, double _sy, double _bz, double _sz) : bx(_bx), sx(_sx), by(_bx), sy(_sy), bz(_bz), sz(_sz) {};
    bool intersect(const Ray& r, double t0, double t1) const;
    void seperate(Node& self, Node& son1, Node& son2, Node& son3, Node& son4, Node& son5, Node& son6, Node& son7, Node& son8);
};


class Surface {
public:
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const = 0;
    virtual bool boxIntersect(Box b);
    virtual void checkPosition(Box& b);
    virtual ~Surface() {};
};

class Triangle : public Surface {  
    SlVector3 a,b,c;
public:
    Triangle(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c) : a(_a), b(_b), c(_c) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    virtual bool boxIntersect(Box b);
    virtual void checkPosition(Box& b);
};

class TrianglePatch : public Triangle {  
    SlVector3 n1, n2, n3;
public: 
    TrianglePatch(const SlVector3 &_a, const SlVector3 &_b, const SlVector3 &_c,
    const SlVector3 &_n1, const SlVector3 &_n2, const SlVector3 &_n3) 
        : Triangle(_a,_b,_c), n1(_n1), n2(_n2), n3(_n3) {};
    virtual bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    //virtual bool boxIntersect(Box b);
    //virtual void checkPosition(Box& b);
};

class Sphere : public Surface {
    SlVector3 c;
    double rad;
public:
    Sphere(const SlVector3 &_c, double _r) : c(_c), rad(_r) {};
    bool intersect(const Ray &r, double t0, double t1, HitRecord &hr) const;
    bool boxIntersect(Box b);
    void checkPosition(Box& b);
};



struct Node {
    Box* self = NULL;
    Node* b000 = NULL;
    Node* b001 = NULL;
    Node* b010 = NULL;
    Node* b011 = NULL;
    Node* b100 = NULL;
    Node* b101 = NULL;
    Node* b110 = NULL;
    Node* b111 = NULL;
};

class Tracer {
    SlVector3 bcolor/*背景颜色*/, eye/*眼睛位置*/, at/*看的方向*/, up/*抬头高度*/;
    double angle, hither;
    unsigned int res[2];
    std::vector<std::pair<Surface *, Fill> > surfaces;
    std::vector<Light> lights;
    double shadowbias;
  
    SlVector3 *im;
public:
    Tracer(const std::string &fname);
    ~Tracer();
    void traceImage();  
    SlVector3 trace(const Ray &ray, double t0, double t1, double c) const;
    SlVector3 shade(const HitRecord &hr, double c) const;
    void writeImage(const std::string &fname);

    bool color;
    int samples;
    double aperture;
    int maxraydepth;
};
	
#endif
