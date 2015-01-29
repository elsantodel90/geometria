#include <iostream>
#include <cmath>

using namespace std;

// INICIO DE CODIGO DE GEOMETRIA

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
    pto proyectar(const pto &p) const { return ((p*d) / d.normaSqr()) * d; } // Proyecta sobre el nuestro
    pto reflejar(const pto &p) const { return p - 2 * normal().proyectar(p); } // Refleja al otro vector usandonos como eje de simetria
    pto simetrico(const pto &p) const { return 2 * (*this) - p; } // Simetrico de p con respecto a o
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
    linea(pto p, pto dd) : s(p), d(dd) {}
    static linea por(const pto &a, const pto &b) { return linea(a, b-a); } // linea::por(a,b)
    linea normalPor(const pto &p) const { return linea(p, d.normal()); }
    linea paralelaPor(const pto &p) const { return linea(p, d); }
    bool contiene(const pto &p) const { return feq(d ^ (p-s),0); }
    floating dist(const pto &p) const { return fabs(distConSigno()); }
    floating distConSigno(const pto &p) const { return (d ^ (p-s)) / d.norma(); }
    pto proyeccion(const pto &p) const { return s + d.proyectar(p-s); }
    pto reflejo(const pto &p) const { return s + d.reflejar(p-s); }
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
    floating dist(const pto &p) const {
        floating h = recta().dist(p), c = hypot(longitud(), h);
        floating d1 = dist(f, p), d2 = dist(t, p);
        if(b<EPSILON || c <= d1 || c <= d2)return min(d1, d2); else return h;
    }

};

enum StatusLineaSegmento {NADA, PARALELO, SE_CORTAN, CONTENIDO};

StatusLineaSegmento status(const linea &l, const segmento &s) {
    linea r = s.recta(); StatusRectas stat = status(l, r);
    if (stat == PARALELAS) return PARALELO;
    if (stat == COINCIDENTES) return CONTENIDO;
    if (s.contiene(interseccion(r, l))) return SE_CORTAN;
    return NADA;
}

tdbl dist(seg a, seg b){
	return (inter(a, b))?0.0:min(min(dist(a.f, b), dist(a.s, b)), min(dist(b.f, a), dist(b.s, a)));
}


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
    pto ortocentro() const { return interseccion(linea::por(a, linea::por(b,c).proyeccion(a)),
                                                 linea::por(b, linea::por(a,c).proyeccion(b))) }
    pto baricentro() const { return (a+b+c) / 3.0; }
    pto incentro() const { return interseccion(bisectriz(a,b,c), bisectriz(b,a, c)); }
    floating inradio() const { return 2 * area() / perimetro(); } // O tambien, linea::por(a,b).dist(incentro())
    floating perimetro() const { return dist(a,b) + dist(b,c) + dist(c,a); }
    floating area() const { return 0.5 * fabs((b-a)^(c-a)); }
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

enum StatusPuntoPoligono {EXTERIOR, FRONTERA, INTERIOR};

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
    pto centroDeMasa() const { // En general, la idea para lograr calcular una integral sobre el poligono es llevarla al borde con teorema de Green.
        escalar A = 0; pto cm(0,0);
        for(int i = 0, j = n-1; i < n; j = i++) {
            A += v[j]^v[i];
            #define green_mix(x,y) (v[j].x*v[j].x + v[i].x*v[j].x + v[i].x*v[i].x) * (v[i].y - v[j].y)
            cm.x += green_mix(x,y); cm.y -= green_mix(y,x);
        }
        cm.x /= (3 * A); cm.y /= (3 * A); return cm; // Notar que hasta aca podian ser todos enteros, y el centro de masa racional.
    }
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
    StatusPuntoPoligono status(const pto &p) const { // Esta escrito para coordenadas enteras, con floats habria que usar comparaciones estandar con EPS
        unsigned i, j, m, M, c = 0;
        for(i = 0, j = v.size()-1; i < v.size(); j = i++){
            if(feq(dist(p,v[i])+dist(p,v[j]),dist(v[i],v[j]))) return FRONTERA; // Frontera debe manejarse aparte de esta forma
            if((v[i].y <= p.y && p.y < v[j].y)||(v[j].y <= p.y && p.y < v[i].y)){
                m = i; M = j; if(v[i].y > v[j].y) swap(m,M);
                c ^= (((v[m]-p)^(v[M]-p)) > 0);
            }
        }
        return c ? INTERIOR : EXTERIOR;
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

// Par de puntos mas cercano.

#define ord(n,a,b) bool n(const pto &p, const pto &q){return ((p.a==q.a)?(p.b<q.b):(p.a<q.a));}
ord(cp_menorEnX,x,y) ord(cp_menorEnY,y,x)

const int CLOSEST_PAIR_MAX_N = 110000;
pto v [CLOSEST_PAIR_MAX_N]; // INPUT
struct ClosestPair {
    pto a,b; escalar d;
    ClosestPair() : d(INF) {}
    ClosestPair(const pto &p1, const pto &p2) : a(p1), b(p2), d(distSqr(p1,p2)) {}
    void check(const pto &p1, const pto &p1) {
        escalar newD = distSqr(p1, p2);
        if (newD < d) { d = newD; a = p1; b = p2; }
    }
}
ClosestPair cpair(int ini, int fin){
  if(fin-ini <= 16) // Tuneable: como mucho impone un costo total N * 16 para los casos base
  {
      ClosestPair ret;
      forsn(i, ini, fin) forsn(j, ini, i) ret.check(vx[i], vx[j]);
      sort(v+ini, v+fin, cp_menorEnY);
      return ret;
  }
  int m = (ini+fin)/2; escalar lineX = v[m].x;
  ClosestPair res = cpair(ini,m); right = cpair(m,fin);
  res.check(right.a, right.b);
  static pto vy[CLOSEST_PAIR_MAX_N];
  merge(v+ini, v+m, v+m, v+fin, vy + ini, cp_menorEnY);
  forsn(i, ini, fin) v[i] = vy[i];
  int T = ini;
  forsn(t, ini, fin)
  if ((vy[t].x - lineX)*(vy[t].x - lineX) < res.d)
  {
      for(int i = T-1; i >= ini && (vy[i].y - vy[t].y)*(vy[i].y - vy[t].y) < res.d; i--)
          res.check(vy[i], vy[t]);
      vy[T++] = vy[t];
  }
  return res;
}
ClosestPair closest_pair(int N){ // Asume que hay al menos dos puntos. Hace aritmetica con enteros si son puntos enteros.
  sort(v, v+N,cp_menorEnX);
  forsn(i, 1, N) if(v[i].x==v[i-1].x && v[i].y==v[i-1].y) return 0;
  return cpair(0,N);
}


// MINIMUM BOUNDING CIRCLE. Tiempo lineal esperado.

circulo mdisk_deIni(vector<pto> v){ //v.size()<=3
  switch(v.size()) {
    case 0: return circulo(pto(0,0), -1);
    case 1: return circulo(v[0], 0);
    case 2: return circulo(0.5 * (v[0] + v[1]) , 0.5 * dist(v[0], v[1])); break;
    default:
        triangulo tri(v[0], v[1], v[2]);
        if(feq(tri.area(), 0)) {assert(false); swap(v[1],v[2]); v.pop_back(); return mdisk_deIni(v);} // Nunca deberia pasar esto... PPP lo tenia asi.
        pto c = tri.circuncentro(); return circulo(c, dist(v[0],c));
  }
}
circulo mdisk_minDisc(const pto *ini,const pto *fin, vector<pto>& pIni){
  const pto *ivp;
  int i,cantP=pIni.size();
  for(ivp=ini,i=0;i+cantP<2 && ivp!=fin;ivp++,i++) pIni.push_back(*ivp);
  circulo r = mdisk_deIni(pIni);
  for(;i>0;i--) pIni.pop_back();
  for(;ivp!=fin;ivp++) if (dist(*ivp, r.c) > r.r){
    pIni.push_back(*ivp);
    if (cantP<2) r=mdisk_minDisc(ini,ivp,pIni);
    else r=mdisk_deIni(pIni);
    pIni.pop_back();
  }
  return r;
}
circulo mdisk_minDisc(vector<pto> ps){ // ESTA ES LA QUE SE USA. Asume puntos distintos. Tiempo esperado lineal.
  random_shuffle(ps.begin(),ps.end()); vector<pto> e;
  return mdisk_minDisc(ps.data(),ps.data()+ps.size(),e);
};


// FIN DE CODIGO DE GEOMETRIA





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



// ********* FIN DE GEOMETRIA_PPP

int main()
{
    return 0;
}
