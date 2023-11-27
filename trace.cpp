#define _USE_MATH_DEFINES
#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "time.h"
//#include <getopt.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#define MAX DBL_MAX
#define MIN DBL_MIN
#endif

static Node mainNode;
static Box mainBox(MIN, MAX, MIN, MAX, MIN, MAX);
static std::vector<Node> firstNodes(8);
static std::vector<std::vector<Node>> secondNodes(8, std::vector<Node>(8));
static std::vector<Box> firstLay(8);
static std::vector<std::vector<Box>> secondLay(8, std::vector<Box>(8));

// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3& a, const SlVector3& b, const SlVector3& c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
        c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) { return x * x; }

bool Surface::boxIntersect(Box b) { return false; }
void Surface::checkPosition(Box& b) {}

bool Triangle::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {

    // Step 1 Ray-triangle test***************************************************************************
    SlVector3 nor = cross(this->b - this->a, this->c - this->a);
    //SlVector3 nor = cross(this->c - this->a, this->b - this->a);
    normalize(nor);
    double x = dot(r.d, nor);
    if (x != 0) {
        double t = (dot(nor, this->a) - dot(nor, r.e)) / (x);
        SlVector3 p = r.e + t * r.d;
        double cross1 = dot(cross(this->b - this->a, p - this->a), cross(this->b - this->a, this->c - this->a));
        double cross2 = dot(cross(this->c - this->a, p - this->a), cross(this->c - this->a, this->b - this->a));
        double cross3 = dot(cross(this->c - this->b, p - this->b), cross(this->c - this->b, this->a - this->b));
        if (dot(cross(this->b - this->a, p - this->a), cross(this->b - this->a, this->c - this->a)) > 0 && dot(cross(this->c - this->a, p - this->a), cross(this->c - this->a, this->b - this->a)) > 0 && dot(cross(this->c - this->b, p - this->b), cross(this->c - this->b, this->a - this->b)) > 0) {
            if (t > t0 && t < t1) {
                //设置HitRecord
                hr.t = t;
                //hr.beta = beta;
                //hr.gamma = gemma;
                //hr.alpha = 1 - beta - gemma;
                hr.p = p;
                hr.n = nor;
                hr.e = 0.0 - r.d;
                return true;
            }
        }
    }
    return false;
}

bool Triangle::boxIntersect(Box b) {
    double maxx = std::max(this->a[0], std::max(this->b[0], this->c[0]));
    double minx = std::min(this->a[0], std::min(this->b[0], this->c[0]));
    double maxy = std::max(this->a[1], std::max(this->b[1], this->c[1]));
    double miny = std::min(this->a[1], std::min(this->b[1], this->c[1]));
    double maxz = std::max(this->a[2], std::max(this->b[2], this->c[2]));
    double minz = std::min(this->a[2], std::min(this->b[2], this->c[2]));
    if (maxx >= b.sx && minx <= b.bx && maxy >= b.sy && miny <= b.by && maxz >= b.sz && minz <= b.bz) {
        return true;
    }
    return false;
}

void Triangle::checkPosition(Box& b) {
    b.sx = std::min(std::min(this->a[0], std::min(this->b[0], this->c[0])), b.sx);
    b.bx = std::max(std::max(this->a[0], std::max(this->b[0], this->c[0])), b.bx);
    b.sy = std::min(std::min(this->a[1], std::min(this->b[1], this->c[1])), b.sy);
    b.by = std::max(std::max(this->a[1], std::max(this->b[1], this->c[1])), b.by);
    b.sz = std::min(std::min(this->a[2], std::min(this->b[2], this->c[2])), b.sz);
    b.bz = std::max(std::max(this->a[2], std::max(this->b[2], this->c[2])), b.bz);
}

bool TrianglePatch::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {
    bool temp = Triangle::intersect(r, t0, t1, hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
        return true;
    }
    return false;
}

/*bool TrianglePatch::boxIntersect(Box b) {
    return Triangle::boxIntersect(b);
}

void TrianglePatch::checkPosition(Box& b) {
    Triangle::checkPosition(b);
}*/

bool Sphere::intersect(const Ray& r, double t0, double t1, HitRecord& hr) const {

    // Step 1 Sphere-triangle test***************************************************************************
    double tmin, tmax;
    double flag = (dot(r.d, r.e) - dot(r.d, this->c)) * (dot(r.d, r.e) - dot(r.d, this->c)) - dot(r.d, r.d) * (dot((r.e - this->c), (r.e - this->c)) - this->rad * this->rad);
    if (flag > 0) {
        tmin = -(sqrt(flag) + 2 * dot(r.d, r.e - this->c) / (2 * dot(r.d, r.d)));
        tmax = (sqrt(flag) - 2 * dot(r.d, r.e - this->c) / (2 * dot(r.d, r.d)));
        if (tmin > t0 && tmin < t1) {
            //设置HitRecord
            hr.t = tmin;
            hr.p = r.e + tmin * r.d;
            hr.n = hr.p - this->c;
            hr.e = 0.0 - r.d;
            normalize(hr.n);
            return true;
        }
        else if (tmin < 0.000001 && tmax > t0 && tmax < t1) {
            hr.t = tmax;
            hr.p = r.e + tmax * r.d;
            hr.n = hr.p - this->c;
            hr.e = 0.0 - r.d;
            normalize(hr.n);
            return true;
        }
    }
    return false;
}

SlVector3 HitRecord::getRefractionD() const {
    double idotn = dot(0.0 - this->e, this->n);
    double k, a;
    double ri = this->f.ior;
    if (idotn) {
        k = 1.0 - ri * ri * (1.0 - idotn * idotn);
        if (k < 0.0) {
            return {0.0, 0.0, 0.0};
        }
        a = ri * idotn - sqrt(k);
    }
    else {
        ri = 1.0 / this->f.ior;
        k = 1.0 - ri * ri * (1.0 - idotn * idotn);
        a = ri * idotn + sqrt(k);
    }
    return (0.0 - this->e) * ri - this->n * a;
}

bool Sphere::boxIntersect(Box b) {
    if (this->c[0] <= (b.bx + this->rad) && this->c[0] >= (b.sx - this->rad) && this->c[1] <= (b.by + this->rad) && this->c[1] >= (b.sy - this->rad) && this->c[2] <= (b.bz + this->rad) && this->c[2] >= (b.sz - this->rad)) {
        return true;
    }
    return false;
}

void Sphere::checkPosition(Box& b) {
    b.sx = std::min(this->c[0] - this->rad, b.sx);
    b.bx = std::max(this->c[0] + this->rad, b.bx);
    b.sy = std::min(this->c[1] - this->rad, b.sy);
    b.by = std::max(this->c[1] + this->rad, b.by);
    b.sz = std::min(this->c[2] - this->rad, b.sz);
    b.bz = std::max(this->c[2] + this->rad, b.bz);
}



bool Box::intersect(const Ray& r, double t0, double t1) const {
    //std::cout << "intersect" << std::endl;
    double ta, tb;
    double tsx, tbx, tsy, tby, tsz, tbz;
    ta = (this->sx - r.e[0]) / r.d[0];
    tb = (this->bx - r.e[0]) / r.d[0];
    tsx = std::min(ta, tb);
    tbx = std::max(ta, tb);
    ta = (this->sy - r.e[1]) / r.d[1];
    tb = (this->by - r.e[1]) / r.d[1];
    tsy = std::min(ta, tb);
    tby = std::max(ta, tb);
    ta = (this->sz - r.e[2]) / r.d[2];
    tb = (this->bz - r.e[2]) / r.d[2];
    tsz = std::min(ta, tb);
    tbz = std::max(ta, tb);
    double tmin = std::max(tsx, std::max(tsy, tsz));
    double tmax = std::min(tbx, std::min(tby, tbz));
    return tmin < tmax && tmax > t0 && tmax < t1;
}

void Box::seperate(Node& node, Node& son000, Node& son001, Node& son010, Node& son011, Node& son100, Node& son101, Node& son110, Node& son111) { //已存在内容，只赋值
    node.self = this;
    double mdx, mdy, mdz;
    mdx = (node.self->bx + node.self->sx) / 2;
    mdy = (node.self->by + node.self->sy) / 2;
    mdz = (node.self->bz + node.self->sz) / 2;
    son000.self->sx = son010.self->sx = son100.self->sx = son110.self->sx = node.self->sx;
    son000.self->sy = son001.self->sy = son100.self->sy = son101.self->sy = node.self->sy;
    son000.self->sz = son001.self->sz = son010.self->sz = son011.self->sz = node.self->sz;
    son001.self->sx = son011.self->sx = son101.self->sx = son111.self->sx = son000.self->bx = son010.self->bx = son100.self->bx = son110.self->bx = mdx;
    son010.self->sy = son011.self->sy = son110.self->sy = son111.self->sy = son000.self->by = son001.self->by = son100.self->by = son101.self->by = mdy;
    son100.self->sz = son101.self->sz = son110.self->sz = son111.self->sz = son000.self->bz = son001.self->bz = son010.self->bz = son011.self->bz = mdz;
    son001.self->bx = son011.self->bx = son101.self->bx = son111.self->bx = node.self->bx;
    son010.self->by = son011.self->by = son110.self->by = son111.self->by = node.self->by;
    son100.self->bz = son101.self->bz = son110.self->bz = son111.self->bz = node.self->bz;
    node.b000 = &son000;
    node.b001 = &son001;
    node.b010 = &son010;
    node.b011 = &son011;
    node.b100 = &son100;
    node.b101 = &son101;
    node.b110 = &son110;
    node.b111 = &son111;
}

Tracer::Tracer(const std::string& fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        getline(in, line);
        switch (line[0]) {
        case 'b': {
            std::stringstream ss(line);
            ss >> ch >> bcolor[0] >> bcolor[1] >> bcolor[2];
            break;
        }

        case 'v': {
            getline(in, line);
            std::string junk;
            std::stringstream fromss(line);
            fromss >> junk >> eye[0] >> eye[1] >> eye[2];

            getline(in, line);
            std::stringstream atss(line);
            atss >> junk >> at[0] >> at[1] >> at[2];

            getline(in, line);
            std::stringstream upss(line);
            upss >> junk >> up[0] >> up[1] >> up[2];

            getline(in, line);
            std::stringstream angless(line);
            angless >> junk >> angle;

            getline(in, line);
            std::stringstream hitherss(line);
            hitherss >> junk >> hither;

            getline(in, line);
            std::stringstream resolutionss(line);
            resolutionss >> junk >> res[0] >> res[1];
            break;
        }

        case 'p': {
            bool patch = false;
            std::stringstream ssn(line);
            unsigned int nverts;
            if (line[1] == 'p') {
                patch = true;
                ssn >> ch;
            }
            ssn >> ch >> nverts;
            std::vector<SlVector3> vertices;
            std::vector<SlVector3> normals;
            for (unsigned int i = 0; i < nverts; i++) {
                getline(in, line);
                std::stringstream ss(line);
                SlVector3 v, n;
                if (patch) ss >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2];
                else ss >> v[0] >> v[1] >> v[2];
                vertices.push_back(v);
                normals.push_back(n);
            }
            bool makeTriangles = false;
            if (vertices.size() == 3) {
                if (patch) {
                    surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                        normals[0], normals[1], normals[2]), fill));
                }
                else {
                    surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                }
            }
            else if (vertices.size() == 4) {
                SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                if (dot(n0, n1) > 0 && dot(n0, n2) > 0 && dot(n0, n3) > 0) {
                    makeTriangles = true;
                    if (patch) {
                        surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2],
                            normals[0], normals[1], normals[2]), fill));
                        surfaces.push_back(std::pair<Surface*, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3],
                            normals[0], normals[2], normals[3]), fill));
                    }
                    else {
                        surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                        surfaces.push_back(std::pair<Surface*, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                    }
                }
                if (!makeTriangles) {
                    std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                }
            }
            break;
        }

        case 's': {
            std::stringstream ss(line);
            SlVector3 c;
            double r;
            ss >> ch >> c[0] >> c[1] >> c[2] >> r;
            surfaces.push_back(std::pair<Surface*, Fill>(new Sphere(c, r), fill));
            break;
        }

        case 'f': {
            std::stringstream ss(line);
            ss >> ch >> fill.color[0] >> fill.color[1] >> fill.color[2] >> fill.kd >> fill.ks >> fill.shine >> fill.t >> fill.ior;
            break;
        }

        case 'l': {
            std::stringstream ss(line);
            Light l;
            ss >> ch >> l.p[0] >> l.p[1] >> l.p[2];
            if (!ss.eof()) {
                ss >> l.c[0] >> l.c[1] >> l.c[2];
                coloredlights = true;
            }
            lights.push_back(l);
            break;
        }

        default:
            break;
        }
    }
    if (!coloredlights) for (unsigned int i = 0; i < lights.size(); i++) lights[i].c = 1.0 / sqrt(lights.size());
    im = new SlVector3[res[0] * res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete[] im;
    for (unsigned int i = 0; i < surfaces.size(); i++) delete surfaces[i].first;
}

void traverse(Ray r, Node n, double t0, double t1, std::vector<int>& v) {
    if (n.b000 == NULL) {
        if (n.self->intersect(r, t0, t1)) {
            //std::cout << "hitbox" << std::endl;
            v.insert(v.end(), n.self->intersectList.begin(), n.self->intersectList.end());
        }
    }
    else {
        if (n.self->intersect(r, t0, t1)) {
            traverse(r, *n.b000, t0, t1, v);
            traverse(r, *n.b001, t0, t1, v);
            traverse(r, *n.b010, t0, t1, v);
            traverse(r, *n.b011, t0, t1, v);
            traverse(r, *n.b100, t0, t1, v);
            traverse(r, *n.b101, t0, t1, v);
            traverse(r, *n.b110, t0, t1, v);
            traverse(r, *n.b111, t0, t1, v);

        }
    }
}


SlVector3 Tracer::shade(const HitRecord& hr, double c) const {
    //std::cout << "shade" << std::endl;
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;
    SlVector3 n;
    bool inside;
    if (dot(hr.e, hr.n) < 0) {
        n = 0.0 - hr.n;
        inside = 1;
    }
    else {
        n = hr.n;
        inside = 0;
    }

    for (unsigned int i = 0; i < lights.size(); i++) {
        const Light& light = lights[i];
        bool shadow = false;
    // Step 3 Check for shadows here
        SlVector3 dir = light.p - hr.p;
        normalize(dir);
        Ray r = Ray(hr.p, dir);
        std::vector<int> hitList;
        traverse(r, mainNode, shadowbias, mag(light.p - hr.p), hitList);
        //std::cout << hitList.size() << " " << "Objects" << std::endl;
        for (unsigned int a = 0; a < hitList.size(); a++) {
            shadow = this->surfaces[hitList[a]].first->intersect(r, shadowbias, mag(light.p - hr.p), dummy);
            if (shadow) {
                break;
            }
        }

        if (!shadow) {

            // Step 2 do shading here
            double spetacular = 1.0;
            double diffuse = 0.0;
            SlVector3 r = getSym(n, dir);

            diffuse = (dot(dir, n) > 0.0) ? dot(dir, n) : 0.0;
            spetacular = (dot(r, hr.e) > 0.0) ? pow(dot(r, hr.e), hr.f.shine) : 0.0;
            //Original Phong BRDF
            color += hr.f.kd * hr.f.color * light.c * diffuse + hr.f.ks * hr.f.color * light.c * spetacular;
        }

    }

    c++;
    double r0 = (1 - hr.f.ior) / (1 + hr.f.ior);
    double R0 = r0 * r0;
    double ar = dot(hr.e, n);
    if (ar > 0) {
        ar = 1 - ar;
    }
    else {
        ar = 1 + ar;
    }
    double R = R0 + (1 - R0) * ar * ar * ar * ar * ar;
    //std::cout << R << std::endl;
    //double R = 1 - dot(hr.e, hr.n);
    SlVector3 reflection = { 0.0, 0.0, 0.0 };
    SlVector3 refraction = { 0.0, 0.0, 0.0 };
    SlVector3 reflectColor = { 0.0, 0.0, 0.0 };
    SlVector3 refractColor = { 0.0, 0.0, 0.0 };
    // Step 4 Add code for computing reflection color here
    if (R != 0) {
        SlVector3 nextReflectionD = getSym(n, hr.e);
        normalize(nextReflectionD);
        Ray nextReflectR = Ray(hr.p, nextReflectionD);
        reflection = trace(nextReflectR, shadowbias, MAX, c);
        reflectColor = reflection * hr.f.kd * hr.f.color * dot(nextReflectR.d, n) + reflection * hr.f.ks * hr.f.color;
    }
    

    // Step 5 Add code for computing refraction color here
    if (c < 2 || !inside) {
        if (R != 1) {
            SlVector3 nextRefractD = hr.getRefractionD();
            if (mag(nextRefractD) != 0) {
                Ray nextRefractR = Ray(hr.p, nextRefractD);
                refraction = trace(nextRefractR, shadowbias, MAX, c);
                refractColor = refraction * hr.f.color;
            }
        }
    }
    
    //color += R * reflectColor + (1 - R) * refractColor;
    color += reflectColor;
    return color;
}



SlVector3 Tracer::trace(const Ray& r, double t0, double t1, double c) const {
    //std::cout << "trace" << std::endl;
    if (c == this->maxraydepth) {
        return { 0.0, 0.0, 0.0 };
    }
    HitRecord tempHR;
    HitRecord hr;
    hr.t = MAX;
    SlVector3 color(bcolor);
    SlVector3 selfColor = SlVector3(0.0);
    bool hit = false;

    // Step 1 See what a ray hits  

    std::vector<int> hitList;
    traverse(r, mainNode, t0, t1, hitList);
    //std::cout << hitList.size() << " " << "Objects" << std::endl;
    for (unsigned int i = 0; i < hitList.size(); i++) {
        //std::cout << hitList[i] << " ";
        if (this->surfaces[hitList[i]].first->intersect(r, t0, t1, tempHR)) {
            hit = true;
            if (tempHR.t < hr.t) {
                hr = tempHR;
                hr.f = this->surfaces[hitList[i]].second;
            }
        }
    }
    /*for (unsigned int i = 0; i < this->surfaces.size(); i++) {
        if (this->surfaces[i].first->intersect(r, t0, t1, tempHR)) {
            hit = true;
            if (tempHR.t < hr.t) {
                hr = tempHR;
                hr.f = this->surfaces[i].second;
            }
        }
    }*/

    if (hit) {
        selfColor += shade(hr, c);
        return selfColor;
    }
    else if (c != 0) {
        return { 0.0, 0.0, 0.0 };
    }
    return color;
}

void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up, w);
    normalize(u);
    SlVector3 v = cross(w, u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI / 180.0) * (angle / 2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3* pixel = im;


    //accelerate process
    for (unsigned int x = 0; x < this->surfaces.size(); x++) {
        this->surfaces[x].first->checkPosition(mainBox);
    }
    for (int x = 0; x < 8; x++) {
        Box b(0,0,0,0,0,0);
        firstLay[x] = b;
        Node n;
        n.self = &firstLay[x];
        firstNodes[x] = n;
    }
    mainBox.seperate(mainNode, firstNodes[0], firstNodes[1], firstNodes[2], firstNodes[3], firstNodes[4], firstNodes[5], firstNodes[6], firstNodes[7]); //第一层赋好了

   
    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            Box b(0, 0, 0, 0, 0, 0);
            secondLay[x][y] = b;
            Node n;
            n.self = &secondLay[x][y];
            secondNodes[x][y] = n;
        }
        firstLay[x].seperate(firstNodes[x], secondNodes[x][0], secondNodes[x][1], secondNodes[x][2], secondNodes[x][3], secondNodes[x][4], secondNodes[x][5], secondNodes[x][6], secondNodes[x][7]);
        for (int y = 0; y < 8; y++) {
            for (int z = 0; z < surfaces.size(); z++) {
                if (this->surfaces[z].first->boxIntersect(secondLay[x][y])) {
                    secondLay[x][y].intersectList.push_back(z);
                }
            }
            
        }
    }

    /*for (int i = 0; i < secondLay[0][7].intersectList.size(); i++) {
        std::cout << secondLay[0][7].intersectList[i] << " ";
    }*/


    for (unsigned int j = 0; j < res[1]; j++) {
        std::cout << j << std::endl;
        for (unsigned int i = 0; i < res[0]; i++, pixel++) {
            //std::cout << i << std::endl;

            SlVector3 result(0.0, 0.0, 0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r - l) * (i + rx) / res[0];
                double y = b + (t - b) * (j + ry) / res[1];
                SlVector3 dir = -d * w + x * u + y * v;

                Ray r(eye, dir);
                normalize(r.d);

                result += trace(r, hither, MAX, 0.0);
            }
            //std::cout << result[0] << " " << result[1] << " " << result[2] << std::endl;
            (*pixel) = result / (samples * (result + 0.3));
        }
    }
}

void Tracer::writeImage(const std::string& fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out << "P6" << "\n" << res[0] << " " << res[1] << "\n" << 255 << "\n";
    SlVector3* pixel = im;
    char val;
    for (unsigned int i = 0; i < res[0] * res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write(&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write(&val, sizeof(unsigned char));
    }
    out.close();
}


int main(int argc, char* argv[]) {
    clock_t start_time, end_time;
    start_time = clock();   //获取开始执行时间
    //int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 2;
    bool color = false;
    /*while ((c = getopt(argc, argv, "a:s:d:c")) != -1) {
        switch(c) {
            case 'a':
            aperture = atof(optarg);
            break;
            case 's':
            samples = atoi(optarg);
            break;
            case 'c':
            color = true;
            break;
            case 'd':
            maxraydepth = atoi(optarg);
            break;
            default:
            abort();
        }
    }

    if (argc-optind != 2) {
        std::cout<<"usage: trace [opts] input.nff output.ppm"<<std::endl;
        for (int i=0; i < argc; i++) std::cout<<argv[i]<<std::endl;
        exit(0);
    }*/

    Tracer tracer("D:\\Coding Relative\\assignment-04\\InputFiles\\refraction test.nnf");
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.traceImage();
    tracer.writeImage("D:\\Coding Relative\\assignment-04\\InputFiles\\refraction test.ppm");
    end_time = clock();     //获取结束时间
    double Times = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("%f seconds\n", Times);
};
