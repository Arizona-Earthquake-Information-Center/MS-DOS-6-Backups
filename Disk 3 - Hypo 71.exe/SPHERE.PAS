{$A+,B-,D+,E+,F-,I+,L+,N+,O-,R-,S+,V+}
{$M 16384,0,655360}

program sphere;                                          { 89 / 1 / 4 }

{ TURBO - Pascal program for computation of optimum fault plane solutions
  Author : Wolfgang Koehler, CIPE, Potsdam, Telegrafenberg, Tel. 310-234
}
label EP,EI,Start,Ende;
{$R+}
{$i INC1 }                        { Include - file , see INC1.PAS }

const NOM = 200; NOM2 = 400;      { Number of fps to store in first step }
type Menge = set of 1..NPmax;
     rec = record Op,Ops,Os1,Os2: integer; ON : real end;
     arr = array[1..NOM2] of real;
     ptr = array[0..NOM2] of integer;
var  delta,smax :  arr;
     Po : ptr;
     Opt   :  array[1..NOM] of rec;
     Pt1,Pt2 : Point;
     PPf,PPs : array[1..NOM2] of Plane;
     P1,P2,P3,P4,P5,PP1 : plane;
     Pfirst,Psec : array[1..2] of Plane;
     all,left,right,plus,minus,top,bot,M1,M2,V1,V2 : Menge;
     event,NOp,i,j,p,ps,s1,s2,s,index : integer;
     t,t3,chi,Nsig,NV,Nr,Nmaxmin,Nmin,dip,rake : real;
     Ch : char; Stan,Typ : string[10];
     head : string[80];
     Nft : boolean;
     FocPar : Text;

function Number(var S : Menge): real;
var p: integer; n: real;
begin n:= 0;
   for p:=1 to NP do if p in S then n:= n + abs(PA[p].sig);
   Number:= n;
end;

procedure Sort(n: integer; var d : arr; var p : ptr);
var i,j,t : integer; di : real;
label M0;
begin
   for i:=1 to n do p[i]:= i;
   for i:=2 to n do
   begin j:= i-1; di:= d[i];
      while d[p[j]] > di do
      begin
         t:= p[j+1]; p[j+1]:= p[j]; p[j]:= t; j:=j-1;
         if j=0 then goto M0;
      end;
M0: end;
end;

begin event:= 0;
   write('Input-data-file (without # and extension (.dat)) :');
   readln(Stan); assign(FocPar,Stan + '.FPS'); rewrite(FocPar);

Start:
   repeat if event > 9 then goto Ende;
      fina:= Stan + chr(event+48);
      assign(stat,fina + '.dat'); {$I-} reset(stat); {$I+}
      inc(event);
   until IOresult = 0;
   readln(stat,head); writeln(FocPar,event,'. ',head);
   writeln(head); NP:= 0;
   repeat NP:= NP + 1; write(NP:3);
       read_Point(NP)
   until eof(stat);

   NOP:= 2*NP; for i:=1 to NOp do Opt[i].ON:= NP;
   plus:= [];
   for p:=1 to NP do if PA[p].sig > 0 then plus:= plus + [p];
   minus:= [1..NP] - plus;

   for p:= 1 to NP do
   begin
      for ps:= p+1 to NP do
      begin Pol(P3,PC[p],PC[ps]); t3:= ssq(P3);
         all :=[1..NP] - [p,ps]; left:= [];
         for s:=1 to NP do if s in all then
         begin t:= PC[s].x*P3.x + PC[s].y*P3.y;
            if (t < -1) or (t > t3) then left:= left + [s];
            Pol(P4,PC[s],P3); Pol(P5,P4,P3);
            delta[s]:= P5.x*P3.y - P5.y*P3.x;
         end;
         right:= all - left;
         delta[p]:= 1E30; delta[ps]:= 1E30;
         Sort(NP,delta,Po); Po[0]:= 0; Po[NP+1]:= 0;
         V1:= left*minus + right*plus;
         Nsig:= Number(all); Nr := Number(V1);
         Nmin:= NP;

         for i:=1 to NP+1 do if Po[i] in all then
         begin
            if Nr < Nsig/2 then NV:= Nr else NV:= Nsig - Nr;
            if NV < Nmin then
            begin Nmin:= NV; s1:= Po[i-1]; s2:= Po[i] end;
            t := abs(PA[Po[i]].sig);
            if Po[i] in V1 then Nr:=Nr - t else Nr:=Nr + t;
         end;

         if Nmin < Opt[NOp].ON then
         begin i:= NOp-1;
            while (Opt[i].ON > Nmin) do
            begin Opt[i+1]:= Opt[i]; i:= i-1;
               if i=0 then goto EI;
            end;
        EI: with Opt[i+1] do
            begin Op:= p; Ops:= ps; Os1:= s1; Os2:= s2; ON:= Nmin;
            end;
         end;
         with Opt[1] do
         writeln(p:3,ps:3,Nsig:7:1,Nmin:7:1,Op:3,Ops:3,Os1:3,Os2:3,ON:7:1);
  EP: end;   {  ps  }
    end;  {  p  }

    for i:=1 to NOp do with Opt[i] do
    begin Pol(P1,PC[Op],PC[Ops]);
       chi:= ssq(P1); PP1.x:= -P1.x/chi; PP1.y:= -P1.y/chi;
       Pfirst[1]:= P1; Pfirst[2]:= P1;
       if Os1=0 then Psec[1]:= PP1 else Pol(Psec[1],PC[Os1],P1);
       if Os2=0 then Psec[2]:= PP1 else Pol(Psec[2],PC[Os2],P1);
       for j:= 1 to 2 do
       begin index:= 2*(i-1) + j;
         Varianz(Pfirst[j],Psec[j],delta[index],smax[index]);
         PPf[index]:= Pfirst[j]; PPs[index]:= Psec[j];
       end;
    end;
    Sort(2*NOp,delta,Po);
    assign(ERG,fina + '.erg'); rewrite(ERG);
    for i:=1 to 40 do
    begin j:= Po[i];  s:= (j +1) div 2;
       Reconv(PPf[j],Pt1); Reconv(PPs[j],Pt2);
       writeln(ERG,delta[j]:6:3,Opt[s].ON:6:1,Pt1.Azi:8:2,Pt1.theta:8:2,
                   Pt2.Azi:8:2,Pt2.theta:8:2);
    end;
    close(ERG); close(stat);
    j:= Po[1]; P1:= PPf[j]; P2:= PPs[j]; Reconv(P1,Pt1); Reconv(P2,Pt2);
    all:= [1..NP];

    { File "XXXX.FPS" :
      Only for the optimum solution : but for all events
      - parameters of nodal planes 1 & 2 : azimuth, theta in degrees
      - variance "delta" , smax,
      - number of stations (NP), number of nonfitting stations
    }
    writeln(FocPar,Pt1.Azi:7:1,Pt1.theta:7:1,Pt2.Azi:14:1,Pt2.theta:6:1,'':7,
                   delta[j]:6:3,smax[j]:7:3,NP:5,Number(all):6:1);
    PP1.x:= 0; PP1.y:= 0; Pol(P3,P1,P2);
    Nft:= AFkt(P1,P2,P3,PP1) < 0;

    if Nft then typ:= 'normal' else typ:= 'thrust';

    { two possible fault planes : strike, dip and rake
      and type of the fault }

    writeln(FocPar,strike(Pt1,Pt2,Nft,dip,rake):7:1,dip:7:1,rake:7:1,
                   strike(Pt2,Pt1,Nft,dip,rake):7:1,dip:7:1,rake:7:1,typ:10);

    Stress_Axes(P1,P2,P4,P5); Reconv(P4,Pt1); Reconv(P5,Pt2);

    { main stress axes : pressure and tension  azimuth and dip }
    { dip = 90 - theta !!! }

    writeln(FocPar,Pt1.Azi:7:1,90-Pt1.theta:7:1,'':7,
                   Pt2.Azi:7:1,90-Pt2.theta:7:1);

    goto Start;  { repeat for the next event number }

 Ende : close(FocPar);

end.
