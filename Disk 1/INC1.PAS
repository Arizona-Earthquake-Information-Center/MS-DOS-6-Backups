{$R+}                                                   { 89 / 1 / 4 }
type plane = record x,y : real end;
                    { projection coord. on the compl. plane }
     Point = record Name: string[4]; Azi,theta,sig: real end;
                    { description record for input data }

const NPmax = 100;             { max. number of stations }
      List_Station : boolean = true;
var  stat,erg,par : text;
     NP : integer;
     fina : string[15];
     PA   : array[1..NPmax] of Point;
     PC   : array[1..NPmax] of plane;

function ssq(var Z : plane) : real;
begin ssq:= sqr(Z.x) + sqr(Z.y) end;

procedure Conv(Pt : Point; var P : plane);
                              { projection sphere --> c. plane }
var chi,t,delta : real;
begin
   delta:= Pi/180*Pt.Azi; t:= Pi/180*Pt.theta; chi:= cos(t)/sin(t);
   if chi < 1E-6 then chi:= 1E-6;
   P.x:= sin(delta)/chi; P.y:= cos(delta)/chi;
end;

procedure ReConv(P : plane; var Pt : Point);
                              { re-projection c. plane --> sphere }
begin with Pt do
   begin Azi:= arctan(P.x/P.y);
     if P.y < 0 then Azi:= Pi+Azi; if Azi < 0 then Azi:= Azi+2*Pi;
     Azi:= 180/Pi*Azi;
     theta:= 180/Pi*arctan(sqrt(ssq(P)));
   end;
end;

procedure read_Point(NP : integer);
begin with PA[NP] do
   begin readln(stat,Name,Azi,theta,sig);
      if List_Station then writeln(Name:5,Azi:9:2,theta:9:2,sig:5:1);
      Conv(PA[NP],PC[NP]);
   end;
end;

procedure Pol(var PP : plane; P1,P2 : plane);
       { comp. of the pol to a great circle given by P1 and P2 }
var Det : real;
begin Det:= P1.x*P2.y - P2.x*P1.y;
   if abs(Det) < 1E-10 then Det:= 1E-10;
   PP.x:= (P1.y - P2.y)/Det; PP.y:= (P2.x - P1.x)/Det;
end;

function skalar(P1,P2 : plane) : real;
        { distance angle between two points given in c-coord }
var s,t1,t2 : real;
begin
   s:= (1+P1.x*P2.x+P1.y*P2.y)/(sqrt(1+ssq(P1))*sqrt(1+ssq(P2)));
   if abs(s) >= 1.0 then skalar:= s/abs(s) else skalar:= s;
end;

procedure Stress_Axes(P1,P2 : Plane; var A1,A2 : Plane);
var P3 : Plane;
    t,t2,chi,c3,s3,x3,xs,ys,yp1,yp2 : real;
begin Pol(P3,P1,P2);
      chi:= sqrt(ssq(P3)); c3:= P3.x/chi; s3:= P3.y/chi;
      x3:= chi; xs:= -1/chi;
      yp1:= P1.y*c3 - P1.x*s3; yp2:= yp1*yp1;
      t:= sqrt(1 + sqr(xs)); t2:= t*t;
      ys:= (-2*yp1*t2 + t*(t2 + yp2))/(yp2 - t2);
      A1.x:= xs*c3 - ys*s3; A1.y:= ys*c3 + xs*s3;
      ys:= (-2*yp1*t2 - t*(t2 + yp2))/(yp2 - t2);
      A2.x:= xs*c3 - ys*s3; A2.y:= ys*c3 + xs*s3;
end;

function Fdelta(P1,P2,P3,Ptst : Plane) : real;
var P4,P5 : Plane;
    Fdel : real;
begin Pol(P4,PTst,P3); Pol(P5,P4,P3);
   Fdel:= skalar(P5,P1); Fdel:= 2*Fdel*sqrt(1-sqr(Fdel));
   if P5.y < P1.y then Fdel:= -Fdel; Fdelta:= Fdel;
end;

function AFkt(P1,P2,P3,PTst : Plane) : real;
begin
   AFKt:= sqrt(1 - sqr(skalar(PTst,P3))) * Fdelta(P1,P2,P3,PTst);
end;

function arccos(x: real) : real;
var ac : real;
begin if abs(x) < 1E-10 then ac:= Pi/2 else ac := arctan(sqrt(1-x*x)/x);
      if x < 0 then ac:= ac + Pi; arccos:= ac;
end;

function strike(P1,P2 : Point; Nft : boolean; var dip,rake : real) : real;
var ff,str,cl : real;
begin str:= P1.azi + 90;
   if str > 360 then str:= str - 360; strike:= str;
   dip:= P1.theta; ff:= Pi/180;
   cl:= sin(P2.theta*ff) * sin(ff*(P2.azi-P1.azi));
   rake:= - arccos(cl)/ff; if not Nft then rake:= 180 + rake;
end;

procedure Varianz(var P1,P2 : plane; var smin,smax : real);
var PP,P3,P4,P5 : plane;
    p : integer;
    sm,sp,Fphi,delta,Fdel,x,y,AF : real;
begin sp:= 0; sm:= 0;
   Pol(P3,P1,P2);
   for p:=1 to NP do
   begin
      AF:= AFkt(P1,P2,P3,PC[p])*PA[p].sig;
      if AF < 0 then sm:= sm + sqr(AF) else sp:= sp + sqr(AF);
   end;
   if sp < sm then
   begin smin:= sp; smax:= sm; P3:= P1; P1:= P2; P2:= P3 end else
   begin smin:= sm; smax:= sp end;
end;
