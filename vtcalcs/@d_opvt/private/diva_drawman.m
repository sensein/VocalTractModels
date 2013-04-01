function diva_drawman(P1,P2,jawh,AV,AVval)
% DIVA_DRAWMAN  Plots the DIVA character with the vocal tract
%   Most of this code is taken from the older C code, with minor modifications.
%   The function takes as input:
%       P1,P2 - Plotting data returned from doAM
%       jawh  - Jaw height parameter
%       AV    - Voicing flag
%       AVval - coordinates for plotting the sin voicing function
%   The jaw jeight parameter is used to position the chin w.r.t. the head. This
%   allows one to to show a fixed jaw position

% Satrajit Ghosh, SpeechLab, Boston University. (c)2001
% $Header: /DIVA.1/classes/@d_opvt/private/diva_drawman.m 2     10/18/01 2:45p Satra $

% $NoKeywords: $

% Setup globals
global RELEASE

% Check hold state of the axis
if ishold,
    holdstate = 1;
else,
    holdstate = 0;
    cla;
end;

hold on;

%V = [140   450   140   650];
%axis(V+[-50 +50 -50 +50]);

eye_x = 315;

line([80 480],[80 680],'color','w');
rectangle('Position',[80 80 400 600],'Facecolor','w');

%/**** main part of head ****/  
line([280, 150],[425 425],'color','k');
line([150, 150],[425 150],'color','k');
line([150, 410],[150 150],'color','k');
line([410, 410],[150 225],'color','k');

%/* nose */
line([410, 450],[225, 300],'color','k');
line([450, 410],[300, 300],'color','k');
line([410, 410],[300, 313],'color','k');

%/**** hair on top of head ****/
%for i=150:5:410,
for i=150:5:280,
    line([i, (i + 5)],[140 150],'color','k');
end;

%/**** hair on back of head ****/
for i=145:5:395,
    line([150, 140],[i, (i + 5)],'color','k');
end;

%/**** eye ****/
%line([eye_x+5, eye_x+45],[200 200]);
%line([eye_x+45, eye_x+45],[200 230]);
%line([eye_x+45, eye_x+5],[230 230]);
%line([eye_x+5, eye_x+5],[230 200]);

rectangle('Position',[eye_x+5,200,40,30],'Facecolor','w');
rectangle('Position',[eye_x+25, 215, 20, 15]);

rectangle('Position',[eye_x+30, 219, 12, 10],'Facecolor','k');

%/**** hair on top of eye ****/
for i=eye_x:5:eye_x+40,
    line([i, (i + 5)],[190 196],'color','k');
end;

%/**** body ****/
rectangle('Position',[250, 475, 100, 100]);

if AV,
     line([275,325],[525,525],'color','k');
     line(AVval(:,1),AVval(:,2),'color','k');
end;

%/* legs */
line([290, 250],[575 650],'color','k');
line([250, 220],[650 650],'color','k');
line([310, 350],[575 650],'color','k');
line([350, 380],[650 650],'color','k');

%/* arms */
line([200, 250],[475 525],'color','k');
line([400, 350],[475 525],'color','k');
%/* neck */
line([280, 280],[475 425],'color','k');
line([320, 320],[475 425],'color','k');
%/**************/

x_off = 450;
y_off =470;
lsc = 10;
ivtx = P1(1,:);
ivty = P1(2,:);
evtx = P1(3,:);
evty = P1(4,:);

%/********vocal tract***********/
line(x_off-lsc*evtx,y_off-lsc*evty,'color','r');
line(x_off-lsc*ivtx,y_off-lsc*ivty,'color','r');

%/* upper lip */
lip_x = evtx(length(evtx));
lip_y = evty(length(evty));
line([x_off-lsc*lip_x,x_off-lsc*lip_x],[y_off-lsc*lip_y,y_off-lsc*lip_y-6],'color','k');
line([x_off-lsc*lip_x,410],[y_off-lsc*lip_y-6,y_off-lsc*lip_y-6],'color','k');

%/* lower lip */
lip_x = ivtx(length(ivtx));
lip_y = ivty(length(ivty));
line([x_off-lsc*lip_x,x_off-lsc*lip_x],[y_off-lsc*lip_y,y_off-lsc*lip_y+7],'color','k');
line([x_off-lsc*lip_x,410],[y_off-lsc*lip_y+7,y_off-lsc*lip_y+7],'color','k');

chin_y = 458.78 + (413.14 - 458.78)*(jawh + 3.0)/6.0;

%/* jaw */

%/* vertical line from bottom of lower lip to chin */
line([410 410],[y_off-lsc*lip_y+7,chin_y],'color','k');

%/* line joining chin to neck */
line([410 320],[chin_y,425],'color','k');

axis ij;

% restore state
if ~holdstate,
    hold off;
end;
