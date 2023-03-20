{$A+,B-,D+,E+,F-,I+,L+,N+,O-,R-,S+,V+}
{$M 16384,0,655360}

program wulff;                                      { 89/ 1 / 4  }

{  graphic display of computed focal plane solutions         }
{  Author : W. Koehler, AdW, CIPE, Potsdam, Tel.  311/234    }

uses Graph,GraphExt,Crt;

{$i INC1 }

procedure Marker(x,y,marktyp,markcol : integer);
var SaveCol : byte;
begin
   SaveCol:= GetColor; SetColor(markcol);
   case marktyp of
   1 : begin line(x-3,y,x+3,y); line(x,y-3,x,y+3); end;
   2 : circle(x,y,3);
   3 : begin SetFillStyle(1,markcol);
         bar(x-3,y-3,x+3,y+3);
       end;
   end;
   SetColor(SaveCol);
end;

procedure ArcW(xm,ym,radius,phi1,phi2 : float);
begin
  EllipseW(xm, ym, 180.0/pi*phi1, 180.0/pi*phi2, radius,radius);
end;

procedure Graphtext(x,y : integer; s : string);
var dx,dy : integer;
begin
    dx:= TextWidth(s); dy:= TextHeight(s);
    SetFillStyle(0,0);  bar(x,y,x+dx,y+dy);
    OutTextXY(x,y,s);
end;
type Menge = set of 1..NPmax;
var Plus,Indiff : Menge;
    Ax1,Ax2,PS1,PS2,P1,P2,P3,P : plane;
    Pt : Point;
    Vx,TestV,Va,t,dt,cr1,cr2,cr3,x,y,x1,y1,x2,y2,
    str1,str2,dip1,dip2,rake1,rake2 : real;
    i,ip,j,ne,hfkt : integer;
    NormFault,StCode : boolean;
    Ch : char; Name : string[77];
    Erg_arr : array[1..100] of record qs,qn : real; Pt1,Pt2 : Point end;

function Zahlstr(x : real; n,m : integer) : string;
var st: string;
begin str(x:n:m,st); ZahlStr:= st end;

procedure wulff_net(P : plane; var x,y : real);
var chi2,del,dd : real;
begin chi2:= 1/ssq(P); dd:= sqrt(chi2 + sqr(chi2)) - chi2;
   x:= dd*P.x; y:= dd*P.y;
end;

procedure Plot_Point(Pt : Point);
{ draws a marker for a single station point in the Wulff net }
var typ,col : integer;
    x,y : real; P : plane;
begin with Pt do
   begin
     if sig > 0 then typ:= 1 else typ:= 2;
     if sig > 0 then col:= 4 else col:= 1;
     if abs(sig) > 1.6 then inc(col,8);
     Conv(Pt,P);
     wulff_net(P,x,y); Marker(x_window(x),y_window(y),typ,col);
     if StCode then OuttextXY(x_window(x)+6,y_window(y)+5,Name);
   end;
end;

procedure Grcircle(P : plane);
{ draws a "great-circle" with pol P  }
var r,sq,phi,dphi,phi1,phi2 : real;
begin
   phi:= arctan(P.x/P.y); if P.y < 0 then phi:= Pi + phi;
   phi:= phi + Pi; if phi > 2*Pi then phi:= phi-2*Pi;
   sq:= sqrt(ssq(P)); r:= sqrt(1+sqr(sq)); dphi:= arctan(1/sq);
   phi1:= phi-dphi; phi2:= phi+dphi;
   if dphi > 0.001 then ArcW(P.x,P.y,r,Pi/2-phi2,Pi/2-phi1) else
     LineW(-cos(phi),sin(phi),cos(phi),-sin(phi));
end;

procedure Clean( Plus : Menge);
{ plots all nonfitting stations }
var p, xp, yp: integer; PP: Plane;
    x,y : real;
begin
   xp:= 10; yp:= 360;
   SetFillStyle(0,0); bar(xp,yp,170,440);
   for p:= 1 to NP do if p in Plus then
   begin Conv(PA[p],PP); Wulff_Net(PP,x,y);
      OutTextXY(xp,yp,PA[p].Name);
      inc(xp,50); if xp > 110 then
      begin xp:= 10; inc(yp,20) end;
      { CircleW(x,y,0.04); }
   end;

end;


function Out(var Plus: Menge; var Nft : boolean; P1,P2 : plane) : real;

{ - computes for all stations, if the polarity reading is in the correct
    quadrant, for a given solution P1 & P2
  - all nonfitting stations build up a failure sum sp, which is the
    function result;
  - the sets "Plus" and "Indiff" contain all nonfitting stations and
    indifferent stations (very nearby nodal planes) resp.
}
var PP,P3,P0 : plane;
    p : integer;
    sm,sp,AF : real;
begin sp:= 0; sm:= 0;
   Pol(P3,P1,P2); Plus:= []; Indiff:= [];
   for p:=1 to NP do
   begin AF:= Afkt(P1,P2,P3,PC[p])* PA[p].sig;
      if AF < 0 then sm:= sm + sqr(AF) else
      begin Plus:= Plus + [p]; sp:= sp + sqr(AF) end;
      if abs(AF) < 0.05 then Indiff:= Indiff + [p];
   end;
   P0.x:= 0; P0.y:= 0; Nft:= AFkt(P1,P2,P3,P0) < 0;
   if sp > sm then begin sp:= sm; Plus:= [1..NP] - Plus end else Nft:= not Nft;
   Plus:= Plus - Indiff;
   Clean(Plus); Out:= sp;
end;

var grdr,grmd,xmax,ymax,x0,y0,dxy,wx,x5,xp1,xp2,ior : integer;
    s : string;
begin
   ClrScr;
   repeat GotoXY(1,1);
      write('Name of input-data-file (without extension ".dat") : ');
      clreol;
      readln(fina);
      assign(stat,fina + '.DAT'); {$I-} reset(stat); {$I-}
      ior:= IOResult;
      if ior <> 0 then write('File not found !'#7);
   until ior = 0;
   readln(stat,Name);

   grdr:= 0;  InitGraph(grdr,grmd,'');
   { '' = directory for "egavga.bgi" and the "*.chr" - files }
   xmax:= GetMaxX; ymax:= GetMaxY;

   graphtext(20,ymax-20,Name);
   x0:= 20; y0:= 20; dxy := ymax-50;
   DefineWindow(1, x0,y0,x0+dxy, y0+dxy,false);
   DefineWindow(2,x0+dxy+5,y0,xmax,ymax,true); x5:= 5;

   Window(1,1,80,1); DirectVideo:= false;
   write('station code in display ? (Y/CR) '); Ch:= readkey;
   StCode:= upcase(Ch) = 'Y';
   List_Station:= false; { no printout of stations while reading file }
   NP:= 0;

   repeat NP:= NP + 1; read_Point(NP); until eof(stat);

   assign(ERG,fina + '.ERG'); {$I-} reset(ERG); {$I+}
   ne:= 0; if IOResult = 0 then
   repeat ne:= ne+1; with Erg_arr[ne] do
      readln(ERG,qs,qn,Pt1.Azi,Pt1.theta,Pt2.Azi,Pt2.theta);
   until eof(ERG);

   SelectWindow(2); SetLinestyle(0,0,1);
   ShowWindow;

   SetTextStyle(0,0,0);
   graphtext(x5,20,'Used Stations :');
   graphtext(x5,40,'Solution Nr   :');
   graphtext(x5,60,'Variance      :');
   graphtext(x5,80,'Poles of Nodal Planes');
   graphtext(x5,90,'---------------------');
   graphtext(x5,110,'Azimuth :');
   graphtext(x5,130,'Delta   :');
   graphtext(x5,155,'Main Stress Axes');
   graphtext(x5,165,'----------------');
   graphtext(x5,185,'Azimuth :');
   graphtext(x5,205,'Dip     :');
   graphtext(x5,230,'Fault Plane Parameters');
   graphtext(x5,240,'----------------------');
   graphtext(x5,260,'Strike  :');
   graphtext(x5,280,'Dip     :');
   graphtext(x5,300,'Rake    :');
   graphtext(x5,330,'Nonfitting Stations');
   graphtext(x5,340,'-------------------');

   j:= 1;

   while (j>0) and (j <= ne) do with Erg_arr[j] do
   begin

      Conv(Pt1,P1); Conv(Pt2,P2); Pol(P3,P1,P2);
      Stress_Axes(P1,P2,Ax1,Ax2);

      SelectWindow(1);
      DefineWorld(-1,-1,1,1);
      SetColor(15);
      SetLineStyle(1,1,1);
      LineW(-1,0,1,0); LineW(0,-1,0,1);
      SetLineStyle(0,0,1); ArcW(0,0,1,0,2*pi);
      SetTextStyle(2,0,4);
      for i:= 1 to NP do Plot_Point(PA[i]);

      Wulff_net(Ax1,x1,y1); Wulff_net(Ax2,x2,y2);
      Marker(x_window(x1),y_window(y1),3,12);
      Marker(x_window(x2),y_window(y2),3,12);

      SetLineStyle(0,0,1); Grcircle(P1); Grcircle(P2);
      SetLineStyle(1,0,1); SetColor(9); Grcircle(P3);

      SelectWindow(2); xp2:= 60; SetTextStyle(0,0,0);

      Va:= Out(Plus,NormFault,P1,P2);

      xp1 := 120;
      graphtext(xp1,20,ZahlStr(NP,7,0));
      graphtext(xp1,40,ZahlStr(j,7,0));
      graphtext(xp1,60,ZahlStr(qs,7,3));

      graphtext(xp2,110,ZahlStr(Pt1.Azi,7,1) + ZahlStr(Pt2.Azi,7,1));
      graphtext(xp2,130,ZahlStr(Pt1.theta,7,1) + ZahlStr(Pt2.theta,7,1));

      ReConv(Ax1,Pt);
      graphtext(xp2,185,ZahlStr(Pt.Azi,7,1));
      graphtext(xp2,205,ZahlStr(90-Pt.theta,7,1));    { Dip = 90ø - theta }
      ReConv(Ax2,Pt);
      graphtext(xp2+56,185,ZahlStr(Pt.Azi,7,1));
      graphtext(xp2+56,205,ZahlStr(90-Pt.theta,7,1));

      str1:= strike(Pt1,Pt2,Normfault,dip1,rake1);
      str2:= strike(Pt2,Pt1,Normfault,dip2,rake2);
      graphtext(xp2,260,ZahlStr(str1,7,1) + ZahlStr(str2,7,1));
      graphtext(xp2,280,ZahlStr(dip1,7,1) + ZahlStr(dip2,7,1));
      graphtext(xp2,300,ZahlStr(rake1,7,1) + ZahlStr(rake2,7,1));


     { writeln(Va:10:3); Fehler-Varianz }

     { das Folgende alles nur, weil die write-Prozedur im Grafik-Modus
       die Hintergrundfarbe nicht richtig setzt !! }

      TextAttr:= 14;
      GotoXY(1,1); write('Which solution to draw [1..',ne,'] (ET => End) : ');
      wx:= WhereX; s:= '';
      repeat
        GotoXY(wx,1); write('':20); GotoXY(wx,1);
        write(s,'_'); GotoXY(WhereX-1,1);
        ch:= readkey;
        case Ch of
          #8: if length(s) > 0 then dec(s[0]);
          '0'..'9': s:= s + ch;
        end;
      until ch = #13;
      val(s,j,i);
      if j <= 0 then Exit;
      SelectWindow(1); ClearViewPort;
   end;

end.