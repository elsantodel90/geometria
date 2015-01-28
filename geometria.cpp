#include <iostream>
#include <cmath>

using namespace std;

// INICIO DE CODIGO DE NOTEBOOK

typedef long long escalar; // Escalar, puede ser de punto flotante, o tambien entero cuando los puntos y operaciones que hagamos nos mantienen en enteros.
typedef long double floating; // Se usa para marcar que estas operaciones si o si utilizan/devuelven numeros de punto flotante, no enteros.

const escalar INF = 0x7FFFFFFFFFFFFFFFLL;
const floating epsilon = 1e-9;
#define feq(a,b) (fabs((a)-(b)) < epsilon) // Se puede cambiar (cuando aplica) por la igualdad exacta si se usan operaciones en enteros.
                                           // Asi como esta deberia producir codigo que funciona igualmente 
                                           // (aunque se pasa a float para hacer la cuenta con epsilon al pedo).
struct pto{
    escalar x,y;
    pto() : x(0), y(0) {}
    pto(escalar xx, escalar yy) : x(xx), y(yy) {}
    pto operator+ (const pto& o) const { return pto(x + o.x, y + o.y); }
    pto operator- (const pto& o) const { return pto(x - o.x, y - o.y); }
    pto operator- () const { return pto(-x, -y); }
    escalar operator* (const pto& o) const { return x*o.x + y*o.y; }
    escalar operator^ (const pto& o) const { return x*o.y - y*o.x; }
    pto operator* (escalar k) const { return pto(k * x, k * y);}
    pto operator/ (floating k) const { return pto(x / k, y / k);}
    escalar normaSqr() const { return (*this)*(*this);} // TODO: Chequear que esto no compile si nos comemos alguno(s) de los asteriscos.
    floating norma() const { return hypot(x,y);} // hypot hace la cuenta de la hipotenusa mejor de lo que uno la hace de la manera obvia.
    pto normal() const { return pto(-y, x); } // El vector girado 90 grados en sentido antihorario, asumiendo sistema de coordenadas usual
    pto unitario() const { return (*this) / norma(); } // Vector unitario en la direccion de p
    pto rotar(const pto &o, const pto &a) const { //rota p con centro en o para que a quede horizontal (hacia la derecha)
        linea l = linea::por(o,a);
        linea n = l.normalPor(o);
        return o + pto(-n.distConSigno(*this) , l.distConSigno(*this));
    }
    pto rotarAngulo(floating angle) const { // En sentido antihorario
        floating c = cos(angle), s = sin(angle);
        return pto(x * c - y * s, x * s + y * c);
    }
    pto rotarAngulo(const pto &o, floating angle) const { return o + (*this - o).rotarAngulo(angle); }
    // usa: define feq, pi //struct pto()
    floating angulo() const { //[0,2pi) , asume que no es el origen
        floating alpha = atan2(y,x);
        if (alpha < 0.0) alpha += 2.0 * M_PI;
        return alpha;
    }
};
pto operator* (escalar a, const pto &p) { return p * k; }
escalar distSqr(const pto &a, const pto &b) { return (a-b).normaSqr();}
floating dist(const pto &a, const pto &b) { return (a-b).norma();}

struct linea{
    pto s, d; //p + lambda*d
    linea() {}
    linea(pto p, pto dd) : p(s), d(dd) {}
    static linea por(const pto &a, const pto &b) const { return linea(a, b-a); } // linea::por(a,b)
    linea normalPor(const pto &p) const { return linea(p, d.normal()); }
    bool contiene(const pto &p) const { return feq(d ^ (p-s),0); }
    floating dist(const pto &p) const { return fabs(distConSigno()); }
    floating distConSigno(const pto &p) const { return (d ^ (p-s)) / d.norma(); }
    pto proyeccion(const pto &p) const { return s + (((p-s)*d) / d.normaSqr()) * d; }
};

enum StatusRectas {PARALELAS, SECANTES, COINCIDENTES};

StatusRectas status(const linea &l1, const linea &l2){
    if(!feq(l1.d^l2.d, 0)) return SECANTES;
    if(feq(l1.d^(l2.s - l1.s), 0)) return COINCIDENTES;
    return PARALELAS;
}

// Asume que las rectas son secantes, sino divide por cero.
pto interseccion(const linea &l1, const linea &l2){
    return l1.s + ( ((l2.s-l1.s)^l2.d) / (l1.d^l2.d) ) * l1.d;
}

struct segmento{
    pto f, t;
    segmento() {}
    segmento(const pto &a, const pto &b) : f(a), t(b) {}
    bool contiene(const pto &p) const { return feq((p-f)^(t-f),0) && (p-f)*(p-t)<=epsilon; }
    floating longitud() const { return dist(f,t); }
    linea recta() const { return linea::por(f,t); }
    linea mediatriz() const { return recta().normalPor(0.5 * (f+t)); }
};

int sgn(escalar a){ if(feq(a,0)) return 0; return (a>0)? 1:-1; }

// Equivale a que status(s1,s2) sea SE_CORTAN o SOLAPADOS. Hace todas las cuentas en enteros, si escalar es entero (notar que no usa epsilon ni feq directamente).
bool hayInterseccion(const segmento &s1, const segmento &s2) {
    #define signo_seg(s1,s2) sgn((s2.f-s1.f)^(s1.t-s1.f))*sgn((s2.t-s1.f)^(s1.t-s1.f))
    int sgn1 = signo_seg(s1,s2), sgn2 = signo_seg(s2,s1);
    if((sgn1 < 0 && sgn2 <= 0)||(sgn1 <= 0 && sgn2 < 0)) return true;
    return sgn1 == 0 && sgn2 == 0 && (s2.contiene(s1.f) || s2.contiene(s1.t) || s1.contiene(s2.f) || s1.contiene(s2.t));
}

// Paralelos = En rectas paralelas distintas
// Alineados = En la misma recta, con interseccion vacia
// Se cortan = Exactamente un punto de interseccion
// Solapados = Infinitos puntos de interseccion (estan en la misma recta y la interseccion es un segmento)
// Nada = Caso restante (estan en rectas secantes pero los segmentos no se intersecan)
enum StatusSegmentos {NADA, PARALELOS, ALINEADOS, SE_CORTAN, SOLAPADOS};

StatusSegmentos status(const segmento &s1, const segmento &s2){
    linea l1 = s1.recta(), l2 = s2.recta();
    StatusRectas status = status(l1,l2);
    if (status == PARALELAS) return PARALELOS;
    if (status == SECANTES)
    {
        pto inter = interseccion(l1,l2);
        return (s1.contiene(inter) && s2.contiene(inter)) ? SE_CORTAN : NADA;
    }
    floating f = (s2.f-s1.f)*(s1.t-s1.f)/norma2(s1.t-s1.f);
    floating t = (s2.t-s1.f)*(s1.t-s1.f)/norma2(s1.t-s1.f);
    if(f > t) swap(f,t);
    floating dd = min(t,(tipo)1.0) - max(f,(tipo)0.0);
    return feq(dd,0)? SE_CORTAN : (dd>0)? SOLAPADOS : ALINEADOS;
}

linea bisectriz(const pto &o,const pto &b,const pto &c) { return linea(o, dist(o,c)*(b-o)+dist(o,b)*(c-o)); }

struct triangulo
{
    pto a,b,c;
    triangulo() {}
    triangulo(const pto &aa,const pto &bb, const pto &cc) : a(aa), b(bb), c(cc) {}
    pto circuncentro() const { return interseccion(segmento(a,b).mediatriz(), segmento(a,c).mediatriz()); }
    floating circunradio() const { return dist(a, circuncentro()); }
    pto ortocentro() const {
        return interseccion(linea::por(a, linea::por(b,c).proyeccion(a)),
                            linea::por(b, linea::por(a,c).proyeccion(b))) }
    pto baricentro() const { return (a+b+c) / 3.0; }
    pto incentro() const { return interseccion(bisectriz(a,b,c), bisectriz(b,a, c)); }
    floating perimetro() const { return dist(a,b) + dist(b,c) + dist(c,a); }
    floating area() const { return 0.5 * fabs((b-a)^(c-a)); }
    floating inradio(pto a, pto b, pto c) const { return 2 * area() / perimetro(); }
};

enum StatusLineaCirculo {EXTERIOR, TANGENTE, SECANTE};

struct circulo{
    pto c; floating r;
    circulo() {}
    circulo(const pto &cc, floating rr) : c(cc), r(rr) {}
    StatusLineaCirculo status(const linea &l) const {
        floating dd = l.dist(c);
        if(feq(dd,c.r)) return TANGENTE;
        return (dd<c.r)? SECANTE:EXTERIOR;
    }
    pair<pto,pto> interseccion(const linea &l) const { // Asume que la recta no es exterior
        pto p = l.proyeccion(c);
        if(status(l)==TANGENTE) return make_pair(p,p);
        floating h = dist(p,c); pto dif = sqrt(r*r-h*h) * l.d.unitario();
        return make_pair(p+dif, p-dif);
    }
    pair<pto,pto> tangencias(const pto &p) const { // Asume que el punto es exterior o tangente al circulo
        floating dd = sqrt(distSqr(p,c.c)-c.r*c.r);
        return interseccion(circulo(p, dd),c);
    }
};

enum StatusCirculos {EXTERIORES, TANGENCIA_EXTERIOR, SECANTES, TANGENCIA_INTERIOR, CONTENIDA, CONCENTRICAS};

int status(const circulo &c1, const circulo &c2){
    floating dd = dist(c1.c, c2.c);
    if(feq(dd,0.0)) return CONCENTRICAS;
    floating R = max(c1.r,c2.r), r = min(c1.r,c2.r);
    if(feq(dd,R+r)) return TANGENCIA_EXTERIOR;
    if(dd > R+r) return EXTERIORES;
    if(feq(dd+r,R)) return TANGENCIA_INTERIOR;
    if((dd + r < R)) return CONTENIDA;
    return SECANTES;
}

pair<pto,pto> interseccion(const circulo &c1, const circulo &c2){ // Asume que son tangentes o secantes
    floating sqdist = distSqr(c1.c, c2.c);
    floating mult = -(c2.r*c2.r-c1.r*c1.r-sqdist)/(2.0*sqdist);
    return c1.interseccion(linea::por(c1.c, c2.c).normalPor(c1.c + mult*(c2.c-c1.c)));
}

struct poligono
{
    vector<pto> v;
    int n;
    poligono() : n(0) {}
    poligono(const vector<pto> &vv) : v(vv), n(vv.size()) {}
    void push(const pto &p) { v.push_back(p); n++; }
    void pop(const pto &p) { v.pop_back(); n--; }
    escalar areaPor2ConSigno() const {
        escalar res = 0;
        for(int i = 0, j = n-1; i < n; j = i++) res += v[j]^v[i];
        return res;
    }    
    escalar areaPor2() const { return abs(areaPor2ConSigno()); }
    bool antihorario() const { return areaPor2ConSigno() > 0; }
    floating ancho() const { // pol convexo y antihorario. Algoritmo lineal: encuentra para cada lado, el vertice mas lejano
        if(n < 3) return 0; // No se banca degenerados (todos alineados) con n > 2 (loop infinito)
        floating res = HUGE_VAL;
        for(int b = 0, a = n-1, cand = 1, far = 0; b < n; a = b++) {
            while(((v[b]-v[a])^(v[cand]-v[far])) >= 0) {far = cand++; cand %= n;}
            res = min(res, linea::por(v[a], v[b]).dist(pol[far]));
        }
        return res;
    }
    escalar diameterSqr() const { //poligono convexo y antihorario. Maxima distancia entre vertices en tiempo lineal.
        if (n < 2) return 0;
        if (n == 2) return distSqr(v[0],v[1]); // No se banca degenerados (todos alineados) con n > 2 (loop infinito)
        escalar res = 0;
        for(int b = 0, a = n-1, cand = 1, far = 0; b < n; a = b++) {
            while(((v[b] - v[a])^(v[cand]-v[far])) >= 0){
                res = max(res, distSqr(v[a], v[far])); 
                far = cand++; cand %= n; 
            }
            res = max(res, distSqr(v[a], v[far]));
        }
        return res;
    }
    bool isconvex() const { // n >= 3
        bool mayor = false, menor = false;
        for(int i = n-1, prev = n-2, next = 0; next < n; prev = i, i = next++) {
            escalar giro = (v[next]-v[i])^(v[prev]-v[i]);
            mayor |= giro > 0; menor |= giro < 0;
        }
        return !(mayor && menor);
    }
};


// CHULL

pto chull_r;
bool menorEnYluegoEnX(const pto &p1, const pto &p2){
  return (p1.y==p2.y)?(p1.x<p2.x):(p1.y<p2.y);
}
bool menorGrahamScan(const pto &p1,const pto &p2){
  tipo ar = (p1-r)^(p2-r);
  return(ar==0)?(dist2(p1,r)<dist2(p2,r)):ar>0;
}
poligono chull(vector<pto> v) { // Devuelve la chull, sin puntos alineados, en sentido antihorario. Asume puntos distintos.
  if(v.size()<3) return pol;
  chull_r = *(min_element(v.begin(),v.end(),menorEnYluegoEnX));
  sort(v.begin(),v.end(), menorGrahamScan);
  int i=0, s; poligono ch;
  while(i<(int)v.size()){
    if(ch.n>1 && ((ch.v[ch.n-1]-ch.v[ch.n-2])^(v[i]-ch.v[ch.n-2]))<=0) ch.pop();
    else ch.push(v[i++]);
  }
  return ch;
}

// FIN DE CODIGO DE NOTEBOOK

\subsection{Convex Hull}
\begin{code}
// usa: define sqr // struct pto(-,^) // funcion dist2
typedef vector<pto> VP;
pto r;
bool men2(const pto &p1, const pto &p2){
  return (p1.y==p2.y)?(p1.x<p2.x):(p1.y<p2.y);
}
bool operator<(const pto &p1,const pto &p2){
  tipo ar = (p1-r)^(p2-r);
  return(ar==0)?(dist2(p1,r)<dist2(p2,r)):ar>0;
  //< clockwise, >counterclockwise
}

VP chull(VP poly){
  if(poly.size()<3) return poly;
  r = *(min_element(poly.begin(),poly.end(),men2));
  sort(poly.begin(),poly.end());
  int i=0, s;
  VP ch;
  while(i<poly.size()){
    s = ch.size();
    if(s>1 && ((ch[s-1]-ch[s-2])^(poly[i]-ch[s-2]))<=0) //<0 conterclockwise, >0 clockwise, sin el = deja puntos colineares
      ch.pop_back();
    else
      ch.push_back(poly[i++]);
  }
  return ch;
}
\end{code}
\subsection{Point in Poly}
\begin{code}
// usa: define feq, sqr // struct pto(-,^) // funcion dist
bool pnpoly(pto p, vector<pto>& v){
    unsigned i, j, m, M, c = 0; // mirame!!!!!!! unsigned  y nada mas???/////////////////
    for(i = 0, j = v.size()-1; i < v.size(); j = i++){
        if(feq(dist(p,v[i])+dist(p,v[j]),dist(v[i],v[j]))) return true;
        if((v[i].y <= p.y && p.y < v[j].y)||(v[j].y <= p.y && p.y < v[i].y)){
            m = i; M = j; if(v[i].y > v[j].y) swap(m,M);
            if(((v[m]-p)^(v[M]-p)) > 0) c^=1;
        }
    }
    return c;
}
\end{code}

// ********* FIN DE GEOMETRIA_GERMAN

// ********* GEOMETRIA_PPP

\section{Geom}
\subsection{Point in Poly}\begin{code}
usa: algorithm, vector
struct pto { tipo x,y; };
bool pnpoly(vector<pto>&v,pto p){
  unsigned i, j, mi, mj, c = 0;
  for(i=0, j = v.size()-1; i< v.size(); j = i++){
    if((v[i].y<=p.y && p.y<v[j].y) ||
       (v[j].y<=p.y && p.y<v[i].y)){
      mi=i,mj=j; if(v[mi].y>v[mj].y)swap(mi,mj);
      if((p.x-v[mi].x) * (v[mj].y-v[mi].y)
       < (p.y-v[mi].y) * (v[mj].x-v[mi].x)) c^=1;
    }
  } return c;
}
\end{code}\subsection{Convex Hull}\begin{code}
usa: algorithm, vector, sqr
tipo pcruz(tipo x1,tipo y1,tipo x2,tipo y2){return x1*y2-x2*y1;}
struct pto {
  tipo x,y;
  tipo n2(pto &p2)const{
    return sqr(x-p2.x)+sqr(y-p2.y);
  }
} r;
tipo area3(pto a, pto b, pto c){
  return pcruz(b.x-a.x,b.y-a.y,c.x-a.x,c.y-a.y);
}
bool men2(const pto &p1, const pto &p2){
  return (p1.y==p2.y)?(p1.x<p2.x):(p1.y<p2.y);
}
bool operator<(const pto &p1,const pto &p2){
  tipo ar = area3(r,p1,p2);
  return(ar==0)?(p1.n2(r)<p2.n2(r)):ar>0;
  //< clockwise, >counterclockwise
}
typedef vector<pto> VP;
VP chull(VP & l){
  VP res = l;  if(l.size()<3) return res;
  r = *(min_element(res.begin(),res.end(),men2));
  sort(res.begin(),res.end());
  tint i=0;VP ch;ch.push_back(res[i++]);ch.push_back(res[i++]);
  while(i<res.size())  // area3 > clockwise, < counterclockwise
    if(ch.size()>1 && area3(ch[ch.size()-2],ch[ch.size()-1],res[i])<=0)
      ch.pop_back();
    else
      ch.push_back(res[i++]);
  return ch;
}
\end{code}
\subsection{Circulo m\'inimo}
\begin{code}
usa: algorithm, cmath, vector, pto (con < e ==)
usa: sqr, dist2(pto,pto), tint
typedef double tipo;
typedef vector<pto> VP;
struct circ { tipo r; pto c; };
#define eq(a,b) (fabs(a-b)<0.00000000000001)
circ deIni(VP v){ //l.size()<=3
  circ r;  sort(v.begin(), v.end()); unique(v.begin(), v.end());
  switch(v.size()) {
    case 0: r.c.x=r.c.y=0; r.r = -1; break;
    case 1: r.c=v[0]; r.r=0; break;
    case 2: r.c.x=(v[0].x+v[1].x)/2.0;
        r.c.y=(v[0].y+v[1].y)/2.0;
        r.r=dist2(v[0], r.c); break;
    default: {
      tipo A = 2.0 * (v[0].x-v[2].x);tipo B = 2.0 * (v[0].y-v[2].y);
      tipo C = 2.0 * (v[1].x-v[2].x);tipo D = 2.0 * (v[1].y-v[2].y);
      tipo R = sqr(v[0].x)-sqr(v[2].x)+sqr(v[0].y)-sqr(v[2].y);
      tipo P = sqr(v[1].x)-sqr(v[2].x)+sqr(v[1].y)-sqr(v[2].y);
      tipo det = D*A-B*C;
      if(eq(det, 0)) {swap(v[1],v[2]); v.pop_back(); return deIni(v);}
      r.c.x = ( D*R-B*P)/det;
      r.c.y = (-C*R+A*P)/det;
      r.r = dist2(v[0],r.c);
    }
  }
  return r;
}
circ minDisc(VP::iterator ini,VP::iterator fin,VP& pIni){
  VP::iterator ivp;
  int i,cantP=pIni.size();
  for(ivp=ini,i=0;i+cantP<2 && ivp!=fin;ivp++,i++) pIni.push_back(*ivp);
  circ r = deIni(pIni);
  for(;i>0;i--) pIni.pop_back();
  for(;ivp!=fin;ivp++) if (dist2(*ivp, r.c) > r.r){
    pIni.push_back(*ivp);
    if (cantP<2) r=minDisc(ini,ivp,pIni);
    else r=deIni(pIni);
    pIni.pop_back();
  }
  return r;
}
circ minDisc(VP ps){ //ESTA ES LA QUE SE USA
  random_shuffle(ps.begin(),ps.end()); VP e;
  circ r = minDisc(ps.begin(),ps.end(),e);
  r.r=sqrt(r.r); return r;
};
\end{code}
\subsection{M\'aximo rect\'angulo entre puntos}
\begin{code}
usa: vector, map, algorithm
struct pto {
  tint x,y ;bool operator<(const pto&p2)const{
    return (x==p2.x)?(y<p2.y):(x<p2.x);
  }
};
bool us[10005];
vector<pto> v;
tint l,w;
tint maxAr(tint x, tint y,tint i){
  tint marea=0;
  tint arr=0,aba=w;
  bool partido = false;
  for(tint j=i;j<(tint)v.size();j++){
    if(x>=v[j].x)continue;
    tint dx = (v[j].x-x);
    if(!partido){
      tint ar = (aba-arr) * dx;marea>?=ar;
    } else {
      tint ar = (aba-y) * dx;marea>?=ar;
      ar = (y-arr) * dx;marea>?=ar;
    }
    if(v[j].y==y)partido=true;
    if(v[j].y< y)arr>?=v[j].y;
    if(v[j].y> y)aba<?=v[j].y;
  }
  return marea;
}
tint masacre(){
  fill(us,us+10002,false);
  pto c;c.x=0;c.y=0;v.push_back(c);c.x=l;c.y=w;v.push_back(c);
  tint marea = 0;
  sort(v.begin(),v.end());
  for(tint i=0;i<(tint)v.size();i++){
    us[v[i].y]=true;
    marea>?=maxAr(v[i].x,v[i].y,i);
  }
  for(tint i=0;i<10002;i++)if(us[i])marea>?=maxAr(0,i,0);
  return marea;
}
\end{code}
\subsection{M\'axima cantidad de puntos alineados}
\begin{code}
usa: algorithm, vector, map, set, forn, forall(typeof)
struct pto {
  tipo x,y;
  bool operator<(const pto &o)const{
    return (x!=o.x)?(x<o.x):(y<o.y);
  }
};
struct lin{
  tipo a,b,c;//ax+by=c
  bool operator<(const lin& l)const{
    return a!=l.a?a<l.a:(b!=l.b?b<l.b:c<l.c);
  }
};
typedef vector<pto> VP;
tint mcd(tint a, tint b){return (b==0)?a:mcd(b, a%b);}
lin linea(tipo x1, tipo y1, tipo x2, tipo y2){
  lin l;
  tint d = mcd(y2-y1, x1-x2);
  l.a = (y2-y1)/d;
  l.b = (x1-x2)/d;
  l.c = x1*l.a + y1*l.b;
  return l;
}
VP v;
typedef map<lin, int> MLI;
MLI cl;
tint maxLin(){
  cl.clear();
  sort(v.begin(), v.end());
  tint m=1, acc=1;
  forn(i, ((tint)v.size())-1){
    acc=(v[i]<v[i+1])?1:(acc+1);
    m>?=acc;
  }
  forall(i, v){
    set<lin> este;
    forall(j, v){
    if(*i<*j||*j<*i)
      este.insert(linea(i->x, i->y, j->x, j->y));
    }
    forall(l, este)cl[*l]++;
  }
  forall(l, cl){
    m>?= l->second;
  }
  return m;
}
\end{code}
%\newpage
\subsection{Centro de masa y area de un pol\'igono}
\begin{code}
usa: vector, forn
struct pto { tint x,y; };
typedef vector<pto> poly;
tint pcruz(tint x1, tint y1, tint x2, tint y2) { return x1*y2-x2*y1; }
tint area3(const pto& p, const pto& p2, const pto& p3) {
  return pcruz(p2.x-p.x, p2.y-p.y, p3.x-p.x, p3.y-p.y);
}
tint areaPor2(const poly& p) {
  tint a = 0; tint l = p.size()-1;
  forn(i,l-1) a += area3(p[i], p[i+1], p[l]);
  return abs(a);
}
pto bariCentroPor3(const pto& p1, const pto& p2, const pto& p3) {
  pto r;
  r.x = p1.x+p2.x+p3.x; r.y = p1.y+p2.y+p3.y;
  return r;
}
struct ptoD { double x,y; };
ptoD centro(const poly& p) {
  tint a = 0; ptoD r; r.x=r.y=0; tint l = p.size()-1;
  forn(i,l-1) {
    tint act = area3(p[i], p[i+1], p[l]);
    pto pact = bariCentroPor3(p[i], p[i+1], p[l]);
    r.x += act * pact.x; r.y += act * pact.y; a += act;
  } r.x /= (3 * a); r.y /= (3 * a); return r;
}
\end{code}
\subsection{Par de puntos mas cercano}
\begin{code}
usa algorithm, vector, tdbl, tint, tipo, INF, forn, cmath
const tint MAX_N = 10010;
struct pto { tipo x,y;} r;
typedef vector<pto> VP;
#define ord(n,a,b) bool n(const pto &p, const pto &q){ \
  return ((p.a==q.a)?(p.b<q.b):(p.a<q.a));}
#define sqr(a) ((a)*(a))
ord(mx,x,y);
ord(my,y,x);
bool vale(const pto &p){return mx(p,r);};
tipo dist(pto a,pto b){return sqr(a.x-b.x)+sqr(a.y-b.y);}
pto vx[MAX_N];
pto vy[MAX_N];
tint N;
tipo cpair(tint ini, tint fin){
  if(fin-ini==1)return INF;
  if(fin-ini==2)return dist(vx[ini], vx[ini+1]);
  vector<pto> y(fin-ini);
  copy(vy+ini, vy+fin, y.begin());
  tint  m = (ini+fin)/2;
  r = vx[m];
  stable_partition(vy+ini, vy+fin, vale);
  tipo d = min(cpair(ini, m), cpair(m, fin));
  vector<pto> w;
  forn(i, y.size())if(sqr(fabs(y[i].x-vx[m].x))<=d)w.push_back(y[i]);
  forn(i,w.size()){
    for(tint j=i+1;(j<(tint)w.size())
      && sqr(fabs(w[i].y-w[j].y))<d;j++){
      d<?=dist(w[i],w[j]);
    }
  }
  return d;
}
tipo closest_pair(){
  sort(vx, vx+N,mx);
  sort(vy, vy+N,my);
  for(tint i=1;i<N;i++){
    if(vx[i].x==vx[i-1].x && vx[i].y==vx[i-1].y)return 0;
  }
  return sqrt(cpair(0,N));
}
\end{code}
\subsection{CCW}
\begin{code}
struct point {tint x, y;};
int ccw(const point &p0, const point &p1, const point &p2){
    tint dx1, dx2, dy1, dy2;
    dx1 = p1.x - p0.x; dy1 = p1.y - p0.y;
    dx2 = p2.x - p0.x; dy2 = p2.y - p0.y;
    if (dx1*dy2 > dy1*dx2) return +1;
    if (dx1*dy2 < dy1*dx2) return -1;
    if ((dx1*dx2 < 0) || (dy1*dy2 < 0)) return -1;
    if ((dx1*dx1+dy1*dy1) < (dx2*dx2+dy2*dy2))return +1;
    return 0;
}
\end{code}
\subsection{Sweep Line}
\begin{code}
struct pto { tint x,y; bool operator<(const pto&p2)const{
  return (y==p2.y)?(x<p2.x):(y<p2.y);
}};
struct slp{ tint x,y,i;bool f; bool operator<(const slp&p2)const{
  if(y!=p2.y)return y<p2.y;
  if(x!=p2.x)return x<p2.x;
  if(f!=p2.f)return f;
  return i<p2.i;
}};
slp p2slp(pto p,tint i){slp q;q.x=p.x;q.y=p.y;q.i=i;return q;}
tint area3(pto a,pto b,pto c){
  return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
}
tint giro(pto a,pto b,pto c){
  tint a3=area3(a,b,c);
  if(a3<0) return -1; if(a3>0)return 1;
  return 0;
}
bool inter(pair<pto,pto> a, pair<pto,pto> b){
  pto p=a.first,q=a.second,r=b.first,s=b.second;
  if(q<p)swap(p,q);if(s<r)swap(r,s);
  if(r<p){swap(p,r);swap(q,s);}
  tint a1=giro(p,q,r),a2=giro(p,q,s);
  if(a1!=0 || a2!=0){
    return (a1!=a2) && (giro(r,s,p)!=giro(r,s,q));
  } else {
    return !(q<r);
  }
}
tint cant_intersec(vector<pair<pto,pto> >&v){
  tint ic=0;
  set<slp> Q; list<tint> T;
  for(tint i=0;i<(tint)v.size();i++){
    slp p1=p2slp(v[i].first,i);slp p2=p2slp(v[i].second,i);
    if(p2<p1)swap(p1,p2);
    p1.f=true;p2.f=false;
    Q.insert(p1);Q.insert(p2);
  }
  while(Q.size()>0){
    slp p = *(Q.begin());Q.erase(p);
    if(p.f){
      for(list<tint>::iterator it=T.begin();it!=T.end();it++)
        if(inter(v[*it],v[p.i]))ic++;
      T.push_back(p.i);
    } else {
      T.erase(find(T.begin(),T.end(),p.i));
    }
  }
  return ic;
}
\end{code}
\subsection{Intersecci'on de segmentos}
\begin{code}
struct pto{tint x,y;};
struct seg{pto f,s;};
tint sgn(tint a){return a;return (a>0)?1:((a<0)?(-1):0);}
tint pc(pto a, pto b, pto o){return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);}
tint pe(pto a, pto b, pto o){return (a.x-o.x)*(b.x-o.x)+(a.y-o.y)*(b.y-o.y);}
bool inter(seg a, seg b){
	tint ka = sgn(pc(a.f, a.s, b.f))*sgn(pc(a.f, a.s, b.s));
	tint kb = sgn(pc(b.f, b.s, a.f))*sgn(pc(b.f, b.s, a.s));
	if(ka<0 && kb<0)return true; //cruza sin tocar
	if(ka==0 && (pe(a.f,a.s,b.f) <= 0 || pe(a.f,a.s,b.s) <= 0))return true; //b tiene un vertice en a
	if(kb==0 && (pe(b.f,b.s,a.f) <= 0 || pe(b.f,b.s,a.s) <= 0))return true; //a tiene un vertice en b
	return false;
}
\end{code}
\subsection{Distancia entre segmentos}
\begin{code}
tdbl dist(pto p, seg s){
	tdbl a = fabs(tdbl(pc(s.f, s.s, p)));
	tdbl b = hypot(s.f.x-s.s.x,s.f.y-s.s.y),h=a/b, c = hypot(b, h);
	tdbl d1 = hypot(s.f.x-p.x,s.f.y-p.y), d2 = hypot(s.s.x-p.x,s.s.y-p.y);
	if(b<1e-10 || c <= d1 || c <= d2)return min(d1, d2); else return h;
}
tdbl dist(seg a, seg b){
	return (inter(a, b))?0.0:min(min(dist(a.f, b), dist(a.s, b)), min(dist(b.f, a), dist(b.s, a)));
}
\end{code}
\subsection{Cuentitas}
\begin{code}
usa: cmath, algorithm, tipo
struct pto{tipo x,y;};
struct lin{tipo a,b,c;};
struct circ{pto c; tipo r;};
#define sqr(a)((a)*(a))
const double PI = (2.0 * acos(0.0));
pto punto(tipo x, tipo y){pto r;r.x=x;r.y=y;return r;}
const pto cero = punto(0,0);
pto suma(pto o, pto s, tipo k){
  return punto(o.x + s.x * k, o.y + s.y * k);
}
pto sim(pto p, pto c){return suma(c, suma(p,c,-1), -1);}
pto ptoMedio(pto a, pto b){return punto((a.x+b.x)/2.0,(a.y+b.y)/2.0);}
tipo pc(pto a, pto b, pto o){
  return (b.y-o.y)*(a.x-o.x)-(a.y-o.y)*(b.x-o.x);
}
tipo pe(pto a, pto b, pto o){
  return (b.x-o.x)*(a.x-o.x)+(b.y-o.y)*(a.y-o.y);
}
#define sqrd(a,b) (sqr(a.x-b.x)+sqr(a.y-b.y))
tipo dist(pto a, pto b){return sqrt(sqrd(a,b));}
//#define feq(a,b) (fabs((a)-(b))<0.000000000001) para interseccion
#define feq(a,b) (fabs((a)-(b))<0.000000001)
tipo zero(tipo t){return feq(t,0.0)?0.0:t;}
bool alin(pto a, pto b, pto c){  return feq(0, pc(a,b,c));}
bool perp(pto a1, pto a2, pto b1, pto b2){
  return feq(0, pe(suma(a1, a2, -1.0), suma(b1, b2, -1.0), cero));
}
bool hayEL(tipo A11, tipo A12, tipo A21, tipo A22){
  return !feq(0.0, A22*A11-A12*A21);
}
pto ecLineal(tipo A11, tipo A12, tipo A21, tipo A22, tipo R1, tipo R2){
  tipo det = A22*A11-A12*A21;
  return punto((A22*R1-A12*R2)/det,(A11*R2-A21*R1)/det);
}
lin linea(pto p1, pto p2){
  lin l;
  l.b = p2.x-p1.x;
  l.a = p1.y-p2.y;
  l.c = p1.x*l.a + p1.y*l.b;
  return l;
}
bool estaPL(pto p, lin l){return feq(p.x * l.a + p.y * l.b, l.c);}
bool estaPS(pto p, pto a, pto b){
  return feq(dist(p,a)+dist(p,b),dist(b,a));
}
lin bisec(pto o, pto a, pto b){
  tipo da = dist(a,o);
  return linea(o, suma(a, suma(b,a,-1.0), da / (da+dist(b,o))));
}
bool paral(lin l1, lin l2){return !hayEL(l1.a, l1.b, l2.a, l2.b);}
bool hayILL(lin l1, lin l2){ //!paralelas || misma
  return !paral(l1,l2)|| !hayEL(l1.a, l1.c, l2.a, l2.c);
}
pto interLL(lin l1, lin l2){//li==l2->pincha
  return ecLineal(l1.a, l1.b, l2.a, l2.b, l1.c, l2.c);
}
bool hayILS(lin l, pto b1, pto b2){
  lin b = linea(b1,b2);
  if(!hayILL(l,b))return false;
  if(estaPL(b1,l))return true;
  return estaPS(interLL(l,b), b1,b2);
}
pto interLS(lin l, pto b1, pto b2){
  return interLL(l, linea(b1, b2));
}
pto interSS(pto a1, pto a2, pto b1, pto b2){
  return interLS(linea(a1, a2), b1, b2);
}
bool hayISS(pto a1, pto a2, pto b1, pto b2){
  if (estaPS(a1,b1,b2)||estaPS(a2,b1,b2)) return true;
  if (estaPS(b1,a1,a2)||estaPS(b2,a1,a2)) return true;
  lin a = linea(a1,a2), b = linea(b1, b2);
  if(!hayILL(a,b))return false;
  if(paral(a,b))return false;
  pto i = interLL(a,b);
  //sale(i);sale(a1);sale(a2);sale(b1);sale(b2);cout << endl;
  return estaPS(i,a1, a2) && estaPS(i,b1,b2);
}
tipo distPL(pto p, lin l){
  return fabs((l.a * p.x + l.b * p.y - l.c)/sqrt(sqr(l.a)+sqr(l.b)));
}
tipo distPS(pto p, pto a1, pto a2){
  tipo aa = sqrd(a1, a2);
  tipo d = distPL(p, linea(a1, a2));
  tipo xx = aa+sqr(d);
  tipo a1a1 = sqrd(a1, p);
  tipo a2a2 = sqrd(a2, p);
  if(max(a1a1, a2a2) > xx){
    return sqrt(min(a1a1, a2a2));
  }else{
    return d;
  }
}
//
pto bariCentro(pto a, pto b, pto c){
  return punto(
    (a.x + b.x + c.x) / 3.0,
    (a.y + b.y + c.y) / 3.0);
}
pto circunCentro(pto a, pto b, pto c){
  tipo A = 2.0 * (a.x-c.x);tipo B = 2.0 * (a.y-c.y);
  tipo C = 2.0 * (b.x-c.x);tipo D = 2.0 * (b.y-c.y);
  tipo R = sqr(a.x)-sqr(c.x)+sqr(a.y)-sqr(c.y);
  tipo P = sqr(b.x)-sqr(c.x)+sqr(b.y)-sqr(c.y);
  return ecLineal(A,B,C,D,R,P);
}
pto ortoCentro(pto a, pto b, pto c){
  pto A = sim(a, ptoMedio(b,c));
  pto B = sim(b, ptoMedio(a,c));
  pto C = sim(c, ptoMedio(b,a));
  return circunCentro(A,B,C);
}
pto inCentro(pto a, pto b, pto c){
  return interLL(bisec(a, b, c), bisec(b, a, c));
}
pto rotar(pto p, pto o, tipo s, tipo c){
  //gira cw un angulo de sin=s, cos=c
  return punto(
    o.x + (p.x - o.x) * c + (p.y - o.y) * s,
    o.y + (p.x - o.x) * -s + (p.y - o.y) * c
  );
}
bool hayEcCuad(tipo a, tipo b, tipo c){//a*x*x+b*x+c=0 tiene sol real?
  if(feq(a,0.0))return false;
  return zero((b*b-4.0*a*c)) >= 0.0;
}
pair<tipo, tipo> ecCuad(tipo a, tipo b, tipo c){//a*x*x+b*x+c=0
  tipo dx = sqrt(zero(b*b-4.0*a*c));
  return make_pair((-b + dx)/(2.0*a),(-b - dx)/(2.0*a));
}
bool adentroCC(circ g, circ c){//c adentro de g sin tocar?
  return g.r > dist(g.c, c.c) + c.r || !feq(g.r, dist(g.c, c.c) + c.r);
}
bool hayICL(circ c, lin l){
  if(feq(0,l.b)){
    swap(l.a, l.b);
    swap(c.c.x, c.c.y);
  }
  if(feq(0,l.b))return false;
  return hayEcCuad(
    sqr(l.a)+sqr(l.b),
    2.0*l.a*l.b*c.c.y-2.0*(sqr(l.b)*c.c.x+l.c*l.a),
    sqr(l.b)*(sqr(c.c.x)+sqr(c.c.y)-sqr(c.r))+sqr(l.c)-2.0*l.c*l.b*c.c.y
  );
}
pair<pto, pto> interCL(circ c, lin l){
  bool sw=false;
  if(sw=feq(0,l.b)){
    swap(l.a, l.b);
    swap(c.c.x, c.c.y);
  }
  pair<tipo, tipo> rc = ecCuad(
    sqr(l.a)+sqr(l.b),
    2.0*l.a*l.b*c.c.y-2.0*(sqr(l.b)*c.c.x+l.c*l.a),
    sqr(l.b)*(sqr(c.c.x)+sqr(c.c.y)-sqr(c.r))+sqr(l.c)-2.0*l.c*l.b*c.c.y
  );
  pair<pto, pto> p(
    punto(rc.first, (l.c - l.a * rc.first) / l.b),
    punto(rc.second, (l.c - l.a * rc.second) / l.b)
  );
  if(sw){
    swap(p.first.x, p.first.y);
    swap(p.second.x, p.second.y);
  }
  return p;
}
bool hayICC(circ c1, circ c2){
  lin l;
  l.a = c1.c.x-c2.c.x;
  l.b = c1.c.y-c2.c.y;
  l.c = (sqr(c2.r)-sqr(c1.r)+sqr(c1.c.x)-sqr(c2.c.x)+sqr(c1.c.y)
    -sqr(c2.c.y))/2.0;
  return hayICL(c1, l);
 
}
pair<pto, pto> interCC(circ c1, circ c2){
  lin l;
  l.a = c1.c.x-c2.c.x;
  l.b = c1.c.y-c2.c.y;
  l.c = (sqr(c2.r)-sqr(c1.r)+sqr(c1.c.x)-sqr(c2.c.x)+sqr(c1.c.y)
    -sqr(c2.c.y))/2.0;
  return interCL(c1, l);
}

// ********* FIN DE GEOMETRIA_PPP

int main()
{
    return 0;
}
