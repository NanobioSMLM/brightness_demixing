%18/11/24

%- The code saves a table with the extension .mat in the same directory that the "DataInit" file.
clear all;
close all;

%Initialization
% To perform the Flux analysis, the .csv file path contening the localizations information should be indicated in "DataInit" without the extension at the end.
DataInit = 'D:\Laurent\Flux code\Code Analysis Flux Article\Tubuline Atto647N Clathrin Atto655 3D DONALD'
Nb_frame = 1;   % Number of frame on which the molecule appears, this is use for saving the final table that keeps only molecule appearing on more than "Nb_frame" frame. This can be keep equal to 1. 
Radius = 100;    %Radius in nm to determine if two detections in consecutive frames are part of the same molecule

%Opening of the .csv file
Data = table2array(readtable(sprintf('%s.csv',DataInit), 'HeaderLines',1));

% Be careful: "frame" shoud corrrepond to the colum where the frame index on which the localization appear, "x" and "y" to the lateral coordinates, 
% "z" to the axial position, "width" to the sigma of the fitted gaussian, "I" to the intensity (number of photons), "precision" to the uncertainty calculated with Cramer-Rao lower bound
frame=Data(:,2); x=Data(:,12); y=Data(:,13); z=Data(:,5);  width=Data(:,14); I=Data(:,15); ratio=Data(:,15); precision=Data(:,19); % DONALD

clear Data;
%% Plot of the intensity histogram
figure;
Hist_ratio=histogram(I,350,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
Occu_ratio=Hist_ratio.Values;
set(gca,'FontSize',18);
axis([0 0.5*max(I) 0 1.05*max(Occu_ratio)])  
xlabel('Intensity (photon)');
ylabel('Frequency');
axis square;
hold on; 
%% Saving intensity histogram
saveas(gcf,[DataInit,'_Histo_Intensity.pdf']);
saveas(gcf,[DataInit,'_Histo_Intensity.png']);
saveas(gcf,[DataInit,'_Histo_Intensity.fig']);
%%
%Data initialization
N_det=length(x);
molec=zeros(N_det,1);               %Number of the molecule
molec_On=zeros(N_det,1);            %Nber frames in which the molecule is detected (including transiant on and off)
molec_x=zeros(N_det,1);
molec_y=zeros(N_det,1);
molec_z=zeros(N_det,1);
molec_frame_On=zeros(N_det,1);      %Frame at which the molecule turn ON
molec_flux=zeros(N_det,1);
molec_ratio=zeros(N_det,1);         %Average spectral ratio of the molecule INCLUDING the first and last transiant frames
molec_flux_first=zeros(N_det,1);    %Flux measured on the transiant first frame when the molecule turn on
molec_flux_last=zeros(N_det,1);     %Flux measured on the transiant last frame when the molecule turn off
molec_flux_small=zeros(N_det,1);    %Average flux calculated on the frames EXCLUDING the first and last transiant (only non zero for Nb_frame=3)
molec_flux_flat=zeros(N_det,1);     %Average flux calculated on the frames INCLUDING the first and last transiant with normalization by the total flux in the frame
molec_I2=zeros(N_det,1);            %intermediate I^2 used to calculate the standard deviation of the flux of a molecule
molec_std_dev=zeros(N_det,1);   
molec_I2_flat=zeros(N_det,1);       %intermediate I^2 used to calculate the standard deviation of the flat flux of a molecule
molec_I2_first=zeros(N_det,1);
molec_I2_last=zeros(N_det,1);
molec_std_dev_flat=zeros(N_det,1);
molec_ratio2=zeros(N_det,1);        %intermediate ratio^2 used to calculate the standard deviation of the spectral ratio of a molecule
molec_std_dev_ratio=zeros(N_det,1);
molec_std_dev_small=zeros(N_det,1); % Standerd deviation for flux without first and last frame
molec_width=zeros(N_det,1);         %Average width of the PSF of a molecule INCLUDING the first and last transiant frames
molec_ratio_max=zeros(N_det,1);     %Spectral ratio in the frame with max flux of the molecule INCLUDING the first and last transiant frames
molec_width_max=zeros(N_det,1);     %Width of the PSF in the frame with max flux  INCLUDING the first and last transiant frames
molec_flux_max=zeros(N_det,1);      %Max flux of the molecule INCLUDING the first and last transiant frames

%Making of the new molecule array
%Initilization of the first frame
cpt_frame=1;

idx_det_ref=find(frame==frame(1));              %idx des detection du 1er frame
frame_cor(idx_det_ref)=cpt_frame;               %adding a reindexing frame colomn to the data  
frame_nb_det(cpt_frame)=sum(idx_det_ref>0); %nber of detected molec in the frame
frame_flux_mean(cpt_frame)=mean(I(idx_det_ref));

molec(idx_det_ref)=idx_det_ref;                 %idx of the molec for the detection in the data (all new for the first frame)
molec_frame_On(idx_det_ref)=cpt_frame;          %frame at which the molec turned on
molec_On(idx_det_ref)=1;                        %Nber of frames the molec is ON
molec_x(idx_det_ref)=x(idx_det_ref);            %pos x
molec_y(idx_det_ref)=y(idx_det_ref);            %pos y
molec_z(idx_det_ref)=z(idx_det_ref);            %pos z
molec_ratio(idx_det_ref)=ratio(idx_det_ref);
molec_flux(idx_det_ref)=I(idx_det_ref);         %intensity
molec_ratio_max(idx_det_ref)=I(idx_det_ref);    %spectral ratio at max intensity
molec_width_max(idx_det_ref)=width(idx_det_ref);   %width at max intensity
molec_flux_max(idx_det_ref)=I(idx_det_ref);
molec_flux_first(idx_det_ref)=I(idx_det_ref);
molec_I2(idx_det_ref)=I(idx_det_ref).^2;        %intensity^2 for the culculation of the std dev
molec_I2_first(idx_det_ref)=I(idx_det_ref).^2;  
molec_width(idx_det_ref)=width(idx_det_ref);    %width of the PSF
molec_flux_flat(idx_det_ref)=I(idx_det_ref)./frame_flux_mean(cpt_frame);         %intensity
molec_I2_flat(idx_det_ref)=(I(idx_det_ref)./frame_flux_mean(cpt_frame)).^2;        %intensity^2 for the calculation of the std dev
molec_ratio2(idx_det_ref)=ratio(idx_det_ref).^2;        %intensity^2 for the culculation of the std dev

ones_det_ref=ones(size(idx_det_ref));
N_molec=max(molec);  %Nber of detected molecules
idx_det_init=max(idx_det_ref)+1; %next detection idx to process

%%
%Cycle the det by frames until finished
while idx_det_init<=N_det, %N_det=length(x), idx_det_init est la liste des indices des détections sur la 1e frame de l'acq
    cpt_frame=cpt_frame+1; %On avance d'une frame
    idx_det=find(frame==frame(idx_det_init)); %new set of idx associated to the new frame
    frame_cor(idx_det)=cpt_frame;             %Re-idxing the frame 
    
    ones_det=ones(size(idx_det));
    
    %measure the distances^2 of all detected molecules between the new
    %frame and the previous one (reference)
    dist2=(x(idx_det_ref)*ones_det'-ones_det_ref*x(idx_det)').^2+(y(idx_det_ref)*ones_det'-ones_det_ref*y(idx_det)').^2;
    [dist2_nbor,idx_nbor]=min(dist2,[],1);  %Takes the idx and the value of the min
    
    %First new detected molecules
    N_new_molec=sum(~(dist2_nbor<Radius^2));                            %Nber of new molec if the min distance is >= Radius
    idx_det_new_molec=idx_det(~(dist2_nbor<Radius^2));                  %idx of the detections in the new frame that are new molecules 
    idx_new_molec=[(N_molec+1):(N_molec+length(idx_det_new_molec))]';   %new idx to be attributed to the new molecules 
    molec(idx_det_new_molec)=idx_new_molec;                             %idx of the new molecule for the detection data   
    frame_nb_det(cpt_frame)=sum(idx_det>0);                             %nber of detection in the new frame
    frame_flux_mean(cpt_frame)=mean(I(idx_det));
    
    %Attribution of the new molecules found in this frame
    molec_frame_On(idx_new_molec)=cpt_frame;
    molec_On(idx_new_molec)=1;  
    molec_x(idx_new_molec)=x(idx_det_new_molec);
    molec_y(idx_new_molec)=y(idx_det_new_molec);
    molec_z(idx_new_molec)=z(idx_det_new_molec);
    molec_flux(idx_new_molec)=I(idx_det_new_molec);
    molec_flux_first(idx_new_molec)=I(idx_det_new_molec);
    molec_I2_first(idx_new_molec)=I(idx_det_new_molec).^2;
    molec_flux_flat(idx_new_molec)=I(idx_det_new_molec)/frame_flux_mean(cpt_frame);        
    molec_I2(idx_new_molec)=I(idx_det_new_molec).^2;
    molec_I2_flat(idx_new_molec)=(I(idx_det_new_molec)./frame_flux_mean(cpt_frame)).^2; 
    molec_width(idx_new_molec)=width(idx_det_new_molec);
    molec_ratio(idx_new_molec)=ratio(idx_det_new_molec);
    molec_ratio2(idx_new_molec)=ratio(idx_det_new_molec).^2;
    molec_ratio_max(idx_new_molec)=I(idx_det_new_molec);    
    molec_width_max(idx_new_molec)=width(idx_det_new_molec);   
    molec_flux_max(idx_new_molec)=I(idx_det_new_molec);
    
    %Second molecules already detected in the previous frame
    idx_old_molec=molec(idx_det_ref(idx_nbor(dist2_nbor<Radius^2)));    %idx of the molec found in the previous frame 
    idx_det_old_molec=idx_det(dist2_nbor<Radius^2);                     %idx of the detections in the new frame that are associated to the molec detected previously 
    molec(idx_det_old_molec)=idx_old_molec;                             %idx of the old molecule for the detection data
    molec_On(idx_old_molec)=molec_On(idx_old_molec)+1;                  %adding one to the nber of ON frames
    molec_ratio(idx_old_molec)=molec_ratio(idx_old_molec)+ratio(idx_det_old_molec);
    molec_x(idx_old_molec)=molec_x(idx_old_molec)+x(idx_det_old_molec);             %Adding the new coord x (to get the average after)
    molec_y(idx_old_molec)=molec_y(idx_old_molec)+y(idx_det_old_molec);             %Adding the new coord y (to get the average after)
    molec_z(idx_old_molec)=molec_z(idx_old_molec)+z(idx_det_old_molec); 
    molec_flux(idx_old_molec)=molec_flux(idx_old_molec)+I(idx_det_old_molec);       %Adding the new I (to get the average after)
    molec_flux_flat(idx_old_molec)=molec_flux_flat(idx_old_molec)+I(idx_det_old_molec)/frame_flux_mean(cpt_frame);        
    molec_I2(idx_old_molec)=molec_I2(idx_old_molec)+I(idx_det_old_molec).^2;        %Adding the new I^2 (to get the average after)
    molec_I2_flat(idx_old_molec)=molec_I2_flat(idx_old_molec)+(I(idx_det_old_molec)./frame_flux_mean(cpt_frame)).^2;        %intensity^2 for the culculation of the std dev
    molec_width(idx_old_molec)=molec_width(idx_old_molec)+width(idx_det_old_molec); %Adding the new PSF width (to get the average after)
    molec_flux_last(idx_old_molec)=I(idx_det_old_molec);
    molec_I2_last(idx_old_molec)=I(idx_det_old_molec).^2;
    molec_ratio2(idx_old_molec)=molec_ratio2(idx_old_molec)+ratio(idx_det_old_molec).^2;        %Adding the new I^2 (to get the average after)
    if (molec_flux_max(idx_old_molec)<I(idx_det_old_molec)),
        molec_flux_max(idx_old_molec)=I(idx_det_old_molec);
        molec_width_max(idx_old_molec)=width(idx_det_old_molec);
        molec_ratio_max(idx_old_molec)=ratio(idx_det_old_molec);
    end;
    
    %new frame -> ref frame
    idx_det_ref=idx_det;
    idx_det_init=max(idx_det_ref)+1;
    ones_det_ref=ones_det;
    N_molec=N_molec+length(idx_det_new_molec); % 
    
    round(100*idx_det_init/N_det)                   %percentage of cycling 
end;
%%
%Calculation of the averages
molec_x(1:N_molec)=molec_x(1:N_molec)./molec_On(1:N_molec); %molec_On(1:N_molec) est le tableau avec pour chaque molécule le nombre de frame sur lequel c'est ON
molec_y(1:N_molec)=molec_y(1:N_molec)./molec_On(1:N_molec);
molec_z(1:N_molec)=molec_z(1:N_molec)./molec_On(1:N_molec);
molec_flux_small(molec_On>2)=(molec_flux(molec_On>2)-molec_flux_first(molec_On>2)-molec_flux_last(molec_On>2))./(molec_On(molec_On>2)-2);
molec_flux(1:N_molec)=molec_flux(1:N_molec)./molec_On(1:N_molec);
molec_flux_flat(1:N_molec)=molec_flux_flat(1:N_molec)./molec_On(1:N_molec);
molec_std_dev(1:N_molec)=(molec_I2(1:N_molec)./molec_On(1:N_molec)-molec_flux(1:N_molec).^2).^0.5;
molec_I2_first(1:N_molec)=molec_flux_first(1:N_molec).^2;
molec_I2_last(1:N_molec)=molec_flux_last(1:N_molec).^2;
molec_I2_small=zeros(N_det,1);
molec_I2_small(molec_On>2)=molec_I2(molec_On>2)-molec_I2_first(molec_On>2)-molec_I2_last(molec_On>2);
molec_std_dev_small(molec_On>2)=(molec_I2_small(molec_On>2)./(molec_On(molec_On>2)-2)-molec_flux_small(molec_On>2).^2).^0.5;
molec_std_dev_flat(1:N_molec)=(molec_I2_flat(1:N_molec)./molec_On(1:N_molec)-molec_flux_flat(1:N_molec).^2).^0.5;
molec_width(1:N_molec)=molec_width(1:N_molec)./molec_On(1:N_molec);
molec_ratio(1:N_molec)=molec_ratio(1:N_molec)./molec_On(1:N_molec);
molec_std_dev_ratio(1:N_molec)=(molec_ratio2(1:N_molec)./molec_On(1:N_molec)-molec_ratio(1:N_molec).^2).^0.5;

%% Preparation for saving (if Nb_frame=1 all the molecule are saved) 
mol_On=molec_On(molec_On>=Nb_frame);
mol_x=molec_x(molec_On>=Nb_frame);
mol_y=molec_y(molec_On>=Nb_frame);
mol_z=molec_z(molec_On>=Nb_frame);
mol_flux=molec_flux(molec_On>=Nb_frame);
mol_flux_small=molec_flux_small(molec_On>=Nb_frame);
mol_std_dev_small=real(molec_std_dev_small(molec_On>=Nb_frame));
mol_std_dev=molec_std_dev(molec_On>=Nb_frame);
mol_width=molec_width(molec_On>=Nb_frame);
mol_ratio=molec_ratio(molec_On>=Nb_frame);
mol_std_dev_ratio=molec_std_dev_ratio(molec_On>=Nb_frame);
mol_flux_max=molec_flux_max(molec_On>=Nb_frame);
mol_width_max=molec_width_max(molec_On>=Nb_frame);
mol_ratio_max=molec_ratio_max(molec_On>=Nb_frame);

%% Calcul taux de rejection
s = 0; % Compte le nombre de localisation qui aparaissent sur 3 frames min
for k = 1:length(mol_On);
    if mol_On(k) > 2; 
        s = s + mol_On(k);
    end
end
rejection_rate =100- 100*s/length(x) %taux de rejection

%% Data saving with same name as initial data file and .mat extension
NameSaveData = sprintf('%s.mat',DataInit);
save(NameSaveData,'mol_x','mol_y','mol_z','mol_flux','mol_std_dev','mol_flux_small','mol_width','mol_ratio','mol_std_dev_ratio','mol_flux_max','mol_width_max','mol_ratio_max','mol_On','rejection_rate','mol_std_dev_small','molec');
%% Histo Flux
figure;  
xlabel('Flux (photon)');
ylabel('Occurences');
axis square;
hold on;
set(gca,'FontSize',18);
histogram(mol_flux_small(mol_flux_small>0),'Facecolor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
%% Sauvegarde Histo Flux
saveas(gcf,[DataInit,'_Histo_Flux_Code Traitement.pdf']);
saveas(gcf,[DataInit,'_Histo_Flux_Code Traitement.png']);
saveas(gcf,[DataInit,'_Histo_Flux_Code Traitement.fig']);