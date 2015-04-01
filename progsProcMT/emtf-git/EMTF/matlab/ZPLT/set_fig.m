%  set_fig:
%  sets up figure, given plotting limits for rho
%  returns figure handle
% Usage:  [hfig] = set_fig(lims)

function [hfig] = set_fig(lims,pltnum)

%  load in size factors appropriate for this system
pltfrm

%  call with two arguments to make multiple plots
if nargin ==2
  size_fac=size_fac2;
else
  size_fac = size_fac1;
end

globinc
%global sym_line_width line_styles symbol_styles rect rect_paper rect_rho rect_ph

sym_line_width = 1.0;
one_dec = 1.2;
xdecs = log10(lims(2)) - log10(lims(1));
one_dec = one_dec*4/xdecs;
ydecs = log10(lims(4)) - log10(lims(3));
paper_width = xdecs*one_dec;
paper_height = ( ydecs + 3 ) * one_dec;
paper_height = min([paper_height,9]);
rect = [0.5,0.5,paper_width,paper_height] * size_fac;
if nargin > 1
   rect(1) = rect(1) + ((0.2+paper_width)*(pltnum-1))*size_fac;
end
rect_paper = [1.,1.,paper_width,paper_height];

rect_rho = [.15,.15+2.3/(ydecs+3),.8,ydecs/(ydecs+3)*.8];
rect_ph = [.15,.15,.8,2/(ydecs+3)*.8];
hfig = figure('Position',rect,'PaperPosition',rect_paper);
marg = 1.25;

%  First half are
%%%   new colors ... better for matlab 5
line_styles = ['g-';'r-';'b-';'m-';'c-';'k-';'y-';...
               'g-';'r-';'b-';'m-';'c-';'k-';'y-'];
symbol_styles = ['go';'rx';'bo';'mx';'co';'kx';'yo';...
                 'g+';'r*';'b+';'m*';'c+';'k*';'y+'];
%%%   old colors ... OK for matlab 4
%line_styles = ['g-';'m-';'c-';'r-';'y-';'b-';'w-';...
%               'g-';'m-';'c-';'r-';'y-';'b-';'w-'];
%symbol_styles = ['go';'mx';'co';'rx';'yo';'bx';'wo';...
%                 'g+';'m*';'c+';'r*';'y+';'b*';'w+'];
Nsym = 7;
