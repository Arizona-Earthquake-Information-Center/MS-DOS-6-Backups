unit GraphExt;
{$N+}
interface

uses Graph;

type float = double;

procedure DefineWorld(xl,yl,xr,yr : float);
function x_window(x : float) : integer;
function y_window(y : float) : integer;
procedure PutPixelW(x,y : float; Col : byte);
procedure LineW(x1,y1,x2,y2 : float);
procedure DrawPolyW(n : integer; var WP);
procedure EllipseW(xm,ym,phi1,phi2,xr,yr : float);
procedure CoordW(dx,dy : float);
procedure FunctW(n : integer; x1,dx : float; var sy);
procedure RectangleW(x1,y1,x2,y2 : float);

procedure DefineWindow(wn,xl,yo,xr,yu : integer; clipping : boolean);
procedure SelectWindow(wn : integer);
procedure AllWindow;
procedure ShowWindow;
procedure OpenWindow(wn : byte);
procedure CloseWindow(wn : byte);
procedure WindowToDisk(FileName : string);
procedure WindowFromDisk(FileName : string; x,y : Word; BitBlt : Byte);

implementation

var World : record xv,yv,  { WK der linken unteren Ecke des Windows }
                   xo,yo,  { "      rechten oberen  "    "    "   }
                   xf,yf : float end;

procedure DefineWorld(xl,yl,xr,yr : float);
var VP : ViewPortType;
begin GetViewSettings(VP);
   with VP,World do
   begin xv:= xl; yv:= yr; xo:= xr; yo:= yl;
      xf:= (x2-x1)/(xr-xl); yf:= (y2-y1)/(yr-yl);
   end;
end;

function x_window(x : float) : integer;
begin x_window := round( World.xf * (x - World.xv) ) end;

function y_window(y : float) : integer;
begin y_window := round( World.yf * (World.yv - y) ) end;

procedure PutPixelW(x,y : float; Col : byte);
begin PutPixel(x_window(x),y_window(y),Col) end;

procedure LineW(x1,y1,x2,y2 : float);
begin Line(x_window(x1),y_window(y1),x_window(x2),y_window(y2)) end;

procedure DrawPolyW(n : integer; var WP);
var vw : array[1..1] of record wx,wy : float end absolute WP;
    x,y,i : integer;
begin for i:=1 to n do with vw[i] do
   begin x:= x_window(wx); y:= y_window(wy);
     if i=1 then MoveTo(x,y) else LineTo(x,y);
   end;
end;

procedure FunctW(n : integer; x1,dx : float; var sy);
var vy : array[1..1] of float absolute sy;
    i,x,y : integer;
begin for i:=1 to n do
   begin x:= x_window(x1 + (i-1)*dx); y:= y_window(vy[i]);
     if i=1 then MoveTo(x,y) else LineTo(x,y);
   end;
end;

procedure EllipseW(xm,ym,phi1,phi2,xr,yr : float);
var t,dt : float;
    i, n : integer;
begin with World do
   begin
     { die interne Ellipse-Procedur von Graph ist zu ungenau !! ::
       Ellipse(x_window(xm),y_window(ym),phi1,phi2,round(xf*xr),round(yf*yr));
     }
     t:= phi1/180*Pi;  dt:= (phi2 - phi1)/180*Pi;
     { exakte Formel fÅr dx < 1 Pixel :: cos(dt/2) < 1 - 1/r  }
     n:= round(sqrt(8*sqrt(abs(xf*xr * yf*yr)))*dt + 2);
     dt:= dt/n;
     MoveTo(x_window(xm + xr*cos(t)), y_window(ym + yr*sin(t)));

     for i:= 1 to n do
     begin t:= t+dt;
        LineTo(x_window(xm + xr*cos(t)), y_window(ym + yr*sin(t)));
     end;
   end;
end;

procedure CoordW(dx,dy : float);
var i,xvz : integer; x,y : float;
    Zahl : string[8];
begin
   with World do
   begin
      for i:= trunc(xv/dx) to trunc(xo/dx) do
      begin if i=0 then SetLineStyle(SolidLn,1,1)
                   else SetLineStyle(DashedLn,1,1);
         x:= i*dx; LineW(x,yv,x,yo);
      end;
      for i:= trunc(yo/dy) to trunc(yv/dy) do
      begin if i=0 then SetLineStyle(SolidLn,1,1)
                   else SetLineStyle(DashedLn,1,1);
         y:= i*dy; LineW(xv,y,xo,y);
      end;
   end;
end;

procedure RectangleW(x1,y1,x2,y2 : float);
begin Rectangle(x_window(x1),y_window(y1),x_window(x2),y_window(y2)) end;

procedure AllWindow;
begin SetViewPort(0,0,GetMaxX,GetMaxY,true) end;

procedure ShowWindow;
var VP : ViewPortType;
begin
   GetViewSettings(VP);
   with VP do Rectangle(0,0,x2-x1,y2-y1);
end;

const WindNr = 10;
var  WindList  : array[1..WindNr] of record
                   VP: ViewPortType; Size : Word; Pt : Pointer
                 end;

procedure DefineWindow(wn,xl,yo,xr,yu : integer; clipping : boolean);
begin
   if (wn < 1) or (wn > WindNr) then Exit;
   with WindList[wn],VP do
   begin x1:= xl; y1:= yo; x2:= xr; y2:= yu;
      clip:= clipping; Size:= 0;
   end;
end;

procedure SelectWindow(wn : integer);
begin if (wn < 1) or (wn > WindNr) then Exit;
   with WindList[wn],VP do SetViewPort(x1,y1,x2,y2,clip)
end;

procedure OpenWindow(wn : byte);
var xb0,xb1 : Word;
begin if (wn < 1) or (wn > WindNr) then Exit;
   with WindList[wn],VP do
   begin
      SelectWindow(wn);
      Size:= ImageSize(x1,y1,x2,y2); GetMem(Pt,Size);
      GetImage(x1,y1,x2,y2,Pt^); ClearViewPort;
   end;
end;

procedure CloseWindow(wn : byte);
begin if (wn < 1) or (wn > WindNr) then Exit;
   with WindList[wn],VP do
   begin PutImage(x1,y1,Pt^,0); FreeMem(Pt,Size) end;
end;

procedure WindowToDisk(FileName : string);
var DiskFi : File;
    Pt : Pointer; y,WS : Word;
    VP : ViewPortType;
begin assign(DiskFi,FileName); rewrite(DiskFi,1);
   GetViewSettings(VP);
   with VP do
   begin WS:= ImageSize(x1,0,x2,0); GetMem(Pt,WS);
      BlockWrite(DiskFi,WS,2);
      for y:= 0 to y2-y1 do
      begin GetImage(0,y,x2-x1,y,Pt^); BlockWrite(DiskFi,Pt^,WS) end;
   end;
   Close(DiskFi); FreeMem(Pt,WS);
end;

procedure WindowFromDisk(FileName : string; x,y : Word; BitBlt : Byte);
var DiskFi : File;
    Pt : Pointer; FS,n,WS : Word;
begin assign(DiskFi,FileName); {$I-} reset(DiskFi,1); {$I-}
   if IOResult <> 0  then
   begin writeln;
      writeln(#7,'Datei ',FileName,' nicht gefunden !'); Exit;
   end;
   BlockRead(DiskFi,WS,2); GetMem(Pt,WS);
   repeat
      BlockRead(DiskFi,Pt^,WS,n); PutImage(x,y,Pt^,BitBlt);
      inc(y);
   until eof(DiskFi);
   Close(DiskFi); FreeMem(Pt,WS);
end;

var i : byte;
begin
   for i:= 1 to WindNr do DefineWindow(i,0,0,GetMaxX,GetMaxY,true);
   DefineWorld(0,0,1,1);
end.
