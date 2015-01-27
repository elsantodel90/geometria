#include <iostream>
#include <cmath>

using namespace std;

// Clase pto: Recomendado cuando hay que ahorrar espacio, no ponerlo en el notebook. Es facil de programar (consiste fundamentalmente de maquinaria de C++)

typedef long long escalar; // Escalar, puede ser de punto flotante, o tambien entero cuando los puntos y operaciones que hagamos nos mantienen en enteros.
typedef long double floating; // Se usa para marcar que estas operaciones si o si utilizan/devuelven numeros de punto flotante, no enteros.

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
    escalar normaSqr() const { return x*x+y*y;}
    floating norma() const { return hypot(x,y);} // hypot hace la cuenta de la hipotenusa mejor de lo que uno la hace de la manera obvia.
    escalar distSqr() const { return (a-b).normaSqr();}
    floating dist() const { return (a-b).norma();}
};
pto operator* (escalar a, const pto &p) { return p * k; }

const floating epsilon = 1e-9;
#define sqr(n) ((n)*(n))

// INICIO DE CODIGO DE NOTEBOOK

#define feq(a,b) (fabs((a)-(b)) < epsilon)


// FIN DE CODIGO DE NOTEBOOK



// ********* GEOMETRIA_GERMAN

\end{code}
\subsection{L\'inea}
\begin{code}
struct line{
    pto s, d; //s + lambda*d
};
// usa: struct line, pto(-)
line Line(pto a, pto b){
    line res;
    res.s = a;
    res.d = b-a;
    return res;
}
// usa: struct pto(), line
pto normal(line l){
    pto n;
    n.x = -l.d.y;
    n.y = l.d.x;
    return n;
}
//usa: define sqr//struct pto(^,-), line // funcion norma
tipo distPL(pto p, line l){//guarda con el signo
    return (l.d^(p-l.s))/norma(l.d);
}
//usa: define sqr//struct pto(*,-,+), line // operator *(tipo,pto) // funcion norma
pto proyPL(pto p, line l){
    return l.s + ((p-l.s)*l.d)/norma2(l.d) * l.d;
}
//usa: define feq// struct pto(-,^), line
int cantInterLL(line l1, line l2){ //infinitas -> 2
    if(!feq(l1.d^l2.d, 0)) return 1;
    if(feq(l1.d^(l2.s - l1.s), 0)) return 2;
    return 0;
}
//usa: struct pto(+,-,^), line  // operator *(tipo,pto)
pto interLL(line l1, line l2){
    return l1.s + ( ((l2.s-l1.s)^l2.d) / (l1.d^l2.d) ) * l1.d;
}
\end{code}
\newpage
\subsection{Segmento}
\begin{code}
struct segment{
    pto f, t;
};
//usa: define feq, epsilon// struct pto(-,^), segment
bool ptoInSeg(pto p, segment s){
    //return feq(dist(p,s.f) + dist(p,s.t), dist(s.f,s.t));
   return feq((p-s.f)^(s.t-s.f),0) && (p-s.f)*(p-s.t)<=epsilon;
}
//usa: define feq
int sgn(tipo a){
    if(feq(a,0)) return 0;
    return (a>0)? 1:-1;
}
//usa: define feq, epsilon// struct pto(-,^), segment // funcion sgn, ptoInSeg
bool hayInterSS(segment s1, segment s2){
    int sgn1 = sgn((s2.f-s1.f)^(s1.t-s1.f))*sgn((s2.t-s1.f)^(s1.t-s1.f));
    int sgn2 = sgn((s1.f-s2.f)^(s2.t-s2.f))*sgn((s1.t-s2.f)^(s2.t-s2.f));
    if((sgn1 < 0 && sgn2 <= 0)||(sgn1 <= 0 && sgn2 < 0)) return true;
    if(sgn1 == 0 && sgn2 == 0)
        return ptoInSeg(s1.f,s2) || ptoInSeg(s1.t,s2) || ptoInSeg(s2.f,s1) || ptoInSeg(s2.t,s1);
    return false;
}
//usa: define feq, epsilon,sqr// struct pto(+,-,^,*), segment, line // funcion cantInterLL,interLL, norma2, ptoInSeg,Line
int cantInterSS(segment s1, segment s2){ //infinitas -> 2
    line l1, l2;
    l1 = Line(s1.f, s1.t);
    l2 = Line(s2.f, s2.t);
    int cant = cantInterLL(l1,l2);
    if(cant == 0) return 0;
    if(cant == 1){
        pto inter = interLL(l1,l2);
        return (ptoInSeg(inter,s1) && ptoInSeg(inter,s2))? 1:0;
    }
    tipo f = (s2.f-s1.f)*(s1.t-s1.f)/norma2(s1.t-s1.f);
    tipo t = (s2.t-s1.f)*(s1.t-s1.f)/norma2(s1.t-s1.f);
    if(f > t) swap(f,t);
    tipo dd = min(t,(tipo)1.0) - max(f,(tipo)0.0);
    return feq(dd,0)? 1 : (dd>0)? 2:0;
}
\end{code}
\newpage
\subsection{Tri\'angulo}
\begin{code}
//usa: difene sqr // struct pto(+,-),line // funcion normal, Line
pto circuncentro(pto a, pto b, pto c){
    line l1, l2;
    l1.s = (1.0/2)*(a+c); l1.d = normal(Line(a,c));
    l2.s = (1.0/2)*(a+b); l2.d = normal(Line(a,b));
    return interLL(l1,l2);
}
//usa: difene sqr // struct pto(-,^) // funcion dist
tipo circunradio(pto a, pto b, pto c){
    tipo sina = abs((b-a)^(c-a))/(dist(a,b)*dist(a,c));
    return dist(b,c)/(2.0*sina);
}
//usa: difene sqr // struct pto(+,-,*,^),line // operator *(tipo,pto) // funcion norma, proyPL, Line,interLL
pto ortocentro(pto a, pto b, pto c){
    line l1 = Line(a, proyPL(a,Line(b,c)));
    line l2 = Line(b, proyPL(b,Line(a,c)));
    return interLL(l1,l2);
}
//usa: difene sqr // struct pto(+,-),line // operator *(tipo,pto) // funcion dist
line bisectriz(pto o, pto b, pto c){
    line l;
    l.s = o;
    l.d = dist(o,c)*(b-o)+dist(o,b)*(c-o);
    return l;
}
//usa: difene sqr // struct pto(^,-) // funcion dist
tipo inradio(pto a, pto b, pto c){
    tipo p = (dist(a,b) + dist(b,c) + dist(c,a));
    return abs(((b-a)^(c-a))/p);
}
\end{code}
\newpage
\subsection{C\'irculo}
\begin{code}
struct circle{
    pto c;
    tipo r;
};
//usa: difene feq, sqr // struct pto(^,-), line, circle // funcion distPL, norma
int cantInterLC(line l, circle c){
    tipo dd = abs(distPL(c.c, l));
    if(feq(dd,c.r)) return 1;
    return (dd<c.r)? 2:0;
}

//usa: difene feq, sqr // struct pto(^,-,+,*), line, circle  // operator *(tipo,pto) // funcion norma, cantInterLC, proyPL,dist
pair<pto,pto> interLC(line l, circle c){
    pto p = proyPL(c.c, l);
    if(cantInterLC(l,c)==1)return make_pair(p,p);
    tipo dd = sqrt(sqr(c.r)-sqr(dist(p,c.c)));
    pto dif = (dd/norma(l.d))*l.d;
    return make_pair(p+dif, p-dif);
}
//usa: difene feq,sqr // struct pto(), circle  // funcion dist
int cantInterCC(circle c1, circle c2){
    tipo dd = dist(c1.c, c2.c);
    tipo R = max(c1.r,c2.r), r = min(c1.r,c2.r);
    if(feq(dd,c1.r + c2.r) || feq(dd+r,R)) return 1;
    return ((dd < c1.r + c2.r) && (dd > R-r))? 2:0;
}
//usa: difene feq,sqr // struct pto(^,-,+,*), circle, line  // operator *(tipo,pto)// funcion norma, cantInterLC, proyPL,dist,dist2,interLC
pair<pto,pto> interCC(circle c1, circle c2){
    tipo sqdist = dist2(c1.c, c2.c);
    tipo mult = -(sqr(c2.r)-sqr(c1.r)-sqdist)/(2.0*sqdist);
    line l;
    l.s = c1.c + mult*(c2.c-c1.c);
    l.d.x = c2.c.y - c1.c.y;
    l.d.y = c1.c.x - c2.c.x;
    return interLC(l,c1);
}
//usa: difene feq,sqr // struct pto(^,-,+,*), circle, line  // operator *(tipo,pto)// funcion norma, cantInterLC, proyPL,dist,dist2,interLC
pair<pto,pto> tangente(pto p, circle c){
    tipo dd = sqrt(dist2(p,c.c)-sqr(c.r));
    circle cc; cc.c = p; cc.r = dd;
    return interCC(cc,c);
}
\end{code}
\newpage
\subsection{Rotaci\'ones y \'Angulos}
\begin{code}
//usa: define sqr //struct pto(^,-,+), line // funcion norma,normal, distPL, Line
pto rotar(pto o, pto a, pto p){ //rota p con centro en o para que a quede horizontal (hacia la derecha)
    pto res;
    line l = Line(o,a);
    line n; n.s = o; n.d = normal(l);
    res.y = distPL(p,l);
    res.x = -1*distPL(p,n);
    return o + res;
}
// usa: struct pto(+)
pto rotarAngulo(pto o, pto p, tipo angle){
    tipo c = cos(angle), s = sin(angle);
    pto res;
    res.x = (p-o).x * c - (p-o).y * s;
    res.y = (p-o).x * s + (p-o).y * c;
    return res + o;
}

#define pi acos(-1)
// usa: define feq, pi //struct pto()
tipo angulo(pto a){//[0,2pi) si a = 0 devuelve 3pi/2
    if(!feq(a.x,0)){
        tipo ang = atan(a.y/a.x);
        if(a.x < 0) ang += pi;
        if(ang < 0) ang += 2*pi;
        return ang;
    }
    return (a.y > 0)? pi/2:(3.0/2)*pi;
}
\end{code}
\newpage
\subsection{Pol\'igonos}
\begin{code}
// usa: struct pto(-,^)
tipo area(vector<pto>& poly){ //poligonos simples
    if(poly.size() < 3) return 0;
    tipo res = 0;
    for(int i = 1; i < poly.size()-1; i++){
        res += (poly[i]-poly[0])^(poly[i+1]-poly[0]);
    }
    return res/2;
}
// usa: define sqr // struct pto(-,^) // funcion dist
tipo diameter(vector<pto>& poly){ //poly.size() > 1, poly convexo y antihorario
    int n = poly.size();
    if(n == 2) return dist(poly[0],poly[1]);
    int i1 = 0, ii1 = n-1, i2 = 1, ii2 = 0;
    tipo maxdist = 0;
    while(i1 < n){
        while(( (poly[i1]-poly[ii1])^(poly[i2]-poly[ii2]) ) >= 0){
            maxdist = max(maxdist, dist(poly[ii1],poly[ii2]));
            ii2 = i2;
            i2 = (i2 < n-1)? i2+1:0;
        }
        maxdist = max(maxdist, dist(poly[ii1],poly[ii2]));
        ii1 = i1++;
    }
    return maxdist;
}

// usa: struct pto(-,^)
bool isconvex(vector<pto>& poly){
    vector<bool>q(2);
    int n = poly.size();
    forn(i,n){
        if(((poly[(i+1)%n]-poly[i])^(poly[(i+n-1)%n]-poly[i])) < 0) q[0] = true;
        else if(((poly[(i+1)%n]-poly[i])^(poly[(i+n-1)%n]-poly[i])) > 0) q[1] = true;
    }
    return !(q[0]&&q[1]);
}
\end{code}
\newpage
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
