clc
clear all
close all

data_f=[];fs=125;f=[];u=1;h=1;total_px=[];total_py=[];total_f=[];
chardata=['simple_signal.mat'];
sig=load (chardata);
ppg_s=sig.total_ppg';
abp_s=sig.total_abp';
[na,nb]=size(ppg_s);
for kk=1:na
    ppg=ppg_s(kk,:);
    abp=abp_s(kk,:);
    [pks_bp,mag_bp]=peakfinder(abp);
    SBP=mean(mag_bp);

    inv=-abp;
    [minind,mag_min]=peakfinder(inv);
    DBP=find_diastol(abp,minind);

    %-----------------ppg feature extraction


    sig_ppg=2*((ppg-min(ppg))/(max(ppg)-min(ppg)))-1;

    lag=15;
   
    s3=sig_ppg(1:end-lag+1);
    s2=sig_ppg(lag:end);

%-------poincare plot
%      figure
%      plot(s2,s3,'-','LineWidth',0.8)
%         hold on
    z=[-1:0.2:1];
    for j=[0:pi/6:pi];
        y=tan(j).*z;
%         plot(z,y,'k-.','LineWidth',0.7)

    end
%            ylim([-1 1])
    %     %-------------------- find crossing points

k=[0:pi/6:pi];
for  kindex=1:length(k)
    px=[];py=[];
    for j=1:length(s2)-1
        % x=s3(1,j);
        %y=s2(1,j);
        m=(s2(1,j+1)-s2(1,j))/(s3(1,j+1)-s3(1,j));
        x_cross=(-m*s3(1,j)+s2(1,j))./(tan(k(kindex))-m);
        y_cross=tan(k(kindex)).*x_cross;

        if  (x_cross<=s3(1,j))&&  (x_cross >= s3(1,j+1))
            px=[px, x_cross];
            py=[py,y_cross];
        end
        if x_cross>=s3(1,j) &&  x_cross<=s3(1,j+1)
            px=[px, x_cross];
            py=[py,y_cross];
        end


    end
%-----------------extract three time serises
        [min_px,min_py,max_px,max_py,min_nx,min_ny,max_nx,max_ny]=extractup_lowband(px,py);
%               plot(max_py,max_px,'g*');
%             plot(min_py,min_px,'r*');
%             plot(max_ny,max_nx,'r*');
%             plot(min_ny,min_nx,'g*');

iner_x=[min_py,max_ny];
iner_y=[min_px,max_nx];
outer_x=[max_py,min_ny];
outer_y=[max_px,min_nx];


X_in(:,kindex)=iner_x;
Y_in(:,kindex)=iner_y;
% [x_in,indexs]=sort(x_in);
% y_in=y_in(indexs);


X_out(:,kindex)=outer_x;
Y_out(:,kindex)=outer_y;
  total_px=[total_px,px];total_py=[total_py,py];
    
end

%------------------------------feature extraction
  [iner_rrvec1,outer_rrvec1,dis_outvec1,dis_invec1,outper1,inper1,poly_vec1,iner_area1,outer_area1,mid_area]=vector_extract(X_out,Y_out,X_in,Y_in)

     f1=numel( total_py);
    f2=mean(total_px);
    f3=mean(total_py);
    f4=std(total_px);
    f5=std(total_py);
    f6=polyarea(s3,s2);
    f7=std(s2);
    f8=std(s3);
    f9=skewness(s2);
    f10= kurtosis(s2);
    f11=skewness(s3);
    f12= kurtosis(s3);
    f13=skewness(total_px);
    f14=skewness(total_py);
    f15= kurtosis(total_px);
    f16= kurtosis(total_py);
    f17=iner_rrvec1';
    f18=outer_rrvec1';
    f19=min(iner_rrvec1);
    f20=min(outer_rrvec1);
    f21=max(iner_rrvec1);
    f22=max(outer_rrvec1);
    f23=mean(iner_rrvec1);
    f24=mean(outer_rrvec1);
    f25=min(dis_outvec1);
    f26=min(dis_invec1');
    f27=max(dis_outvec1);
    f28=max(dis_invec1');
    f29=mean(dis_outvec1);
    f30=mean(dis_invec1');
    f31=std(dis_outvec1);
    f32=std(dis_invec1');
    f33=min(poly_vec1);
    f34=max(poly_vec1);
  f35=mean(poly_vec1);
    f36=std(poly_vec1);
    f37=entropy(poly_vec1);
    f38=entropy(total_px);
    f39=entropy(total_py);
    f40=iner_area1/outer_area1;
      total_px=[];
    total_py=[];
    feature(:,kk)=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13;f14;f15;f16;f17;f18;f19;f20;f21;f22;f23;f24;f25;f26;f27;f28;f29;f30;f31;f32;f33;f34;f35;f36;f37;f38;f39;f40;dis_outvec1';dis_invec1';outper1;inper1;poly_vec1';iner_area1;outer_area1;mid_area1;SBP;DBP];
% hold off

end
total_f=[total_f,feature];




%    save('f_total.mat','total_f')