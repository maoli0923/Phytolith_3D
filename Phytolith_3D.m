%This is the feature extraction for 3D Phytolith that used in the
%manuscript "ELONGATE DENDRITIC phytoliths as markers for cereal
%domestication and cultivation: a 3D morphometric approach" 
% by Rosalie Hermans, Mao Li, William H. Brightly, Timothy J. Gallaher, Wouter Smagghe, Hannah Lee, Leticia Arco, Lara Stas, Perseverence Savieri, Luc Vrydaghs, Karin Nys, Christophe Snoeck, Caroline A. E. Strömberg

%Contact Mao Li (mli@danforthcenter.org) for any question about this code

%For geometric features
file=dir('*.ply');
for idx=1:length(file)
    [F V]=ply_read(file(idx).name,'tri');
    V=V';F=F';
    boundary = select_holes_and_boundary(V,F);%examine whether there are holes, need mesh processing toolbox
    if ~isempty(boundary)
        F=fill_mesh_holes(V,F,boundary,'closed',500);%fill the holes
    end
    [Vol,Area] = mao_mesh2volume(V,F); %measure volume and surface area
    [K CHV]=convhulln(V);%convex hull volume
    Solidity=Vol/CHV;
    V=V-repmat(mean(V),size(V,1),1);
    [COEFF SCORE latent]=pca(V);
    MaxLength=max(SCORE(:,1))-min(SCORE(:,1)); 
    PC1_sd=std(SCORE(:,1));
    PC2_sd=std(SCORE(:,2));
    PC3_sd=std(SCORE(:,3));
    Elongation=sqrt(latent(2)/latent(1));
    Flatness=sqrt(latent(3)/latent(2));
    threshold=min(SCORE(:,1)):(max(SCORE(:,1))-min(SCORE(:,1)))/3:max(SCORE(:,1));
    for s=1:3
        w=find(SCORE(:,1)>=threshold(s)&SCORE(:,1)<threshold(s+1));
        d(s)=max(SCORE(w,2))-min(SCORE(w,2));
        h(s)=max(SCORE(w,3))-min(SCORE(w,3));
    end
    if Elongation<0.5
        MaxWidth=max(d);
        MaxHeight=max(h);
    else
        MaxWidth=max(SCORE(:,2))-min(SCORE(:,2));
        MaxHeight=max(SCORE(:,3))-min(SCORE(:,3));
    end
    Sphericity=6*pi*(Vol).^2/(Area).^3;
    Feature=[Vol CHV Area Solidity Sphericity MaxLength MaxWidth MaxHeight PC1_sd PC2_sd PC3_sd Elongation Flatness];
    save([file(idx).name(1:end-4) '_Feature.mat'],"Feature");
    
    %for core body
    clear U
    U(:,1)=SCORE(:,2);
    U(:,2)=SCORE(:,3);
    U(:,3)=SCORE(:,1);
    if Elongation<0.5
        partition=3;
        areabreaks = linspace(min(U(:,3)),max(U(:,3)),partition+1);
        for s=1:partition
            a=find(U(:,3)>=areabreaks(s)&U(:,3)<=areabreaks(s+1));
            U(a,1:2)=U(a,1:2)-mean(U(a,1:2));
        end
    end
    [theta,rho,z] = cart2pol(U(:,1),U(:,2),U(:,3));
    M=ceil(max(U(:,1)));
    [x,y,z] = pol2cart(theta,max(M-rho,0),z);
    tip1=0.5;tip2=0.5;
    a=find(U(:,3)>max(U(:,3))-tip1*M);
    U2=U(a,:)-repmat([0 0 max(U(:,3))-tip1*M],length(a),1);
    [azimuth,elevation,r] = cart2sph(U2(:,1),U2(:,2),U2(:,3));
    [x2,y2,z2] = sph2cart(azimuth,elevation,max(M-r,0));
    x(a)=x2;y(a)=y2;z(a)=z2+max(U(:,3))-tip1*M;
    a=find(U(:,3)<min(U(:,3))+tip2*M);
    U2=U(a,:)-repmat([0 0 min(U(:,3))+tip2*M],length(a),1);
    [azimuth,elevation,r] = cart2sph(U2(:,1),U2(:,2),U2(:,3));
    [x2,y2,z2] = sph2cart(azimuth,elevation,max(M-r,0));
    x(a)=x2;y(a)=y2;z(a)=z2+min(U(:,3))+tip2*M;
    clear boundary
    K=boundary(x,y,z);
    [IBDV,IBDS] = mao_mesh2volume(V,K);
    Body=IBDV/Vol;
    index=unique(K(:));
    V=V(index,:);
    [K CHV]=convhulln(V);%convex hull volume
    Solidity=IBDV/CHV;
    Sphericity=6*pi*(IBDV).^2/(IBDS).^3;
    V=V-repmat(mean(V),size(V,1),1);
    [COEFF SCORE latent]=pca(V);
    MaxLength=max(SCORE(:,1))-min(SCORE(:,1)); %maximum length
    PC1_sd=std(SCORE(:,1));
    PC2_sd=std(SCORE(:,2));
    PC3_sd=std(SCORE(:,3));
    Elongation=sqrt(latent(2)/latent(1));
    Flatness=sqrt(latent(3)/latent(2));
    threshold=min(SCORE(:,1)):(max(SCORE(:,1))-min(SCORE(:,1)))/3:max(SCORE(:,1));
    for s=1:3
    w=find(SCORE(:,1)>=threshold(s)&SCORE(:,1)<threshold(s+1));
    d(s)=max(SCORE(w,2))-min(SCORE(w,2));
    h(s)=max(SCORE(w,3))-min(SCORE(w,3));
    end
    if Elongation<0.5
    MaxWidth=max(d);
    MaxHeight=max(h);
    else
        MaxWidth=max(SCORE(:,2))-min(SCORE(:,2));
        MaxHeight=max(SCORE(:,3))-min(SCORE(:,3));
    end
    CoreFeature=[IBDV CHV IBDS Body Solidity Sphericity MaxLength MaxWidth MaxHeight PC1_sd PC2_sd PC3_sd Elongation Flatness];
    save([file(idx).name(1:end-4) '_CoreFeature.mat'],"CoreFeature");    

    
%for branching features, need to download javaplex
%(https://appliedtopology.github.io/javaplex/)
    Source=unique(K(:));
    E=[[F(:,1) F(:,2)];[F(:,1) F(:,3)];[F(:,2) F(:,3)]];
    E=[E;[E(:,2) E(:,1)]];
    E=unique(E,'rows');
    a=find(E(:,1)>E(:,2));
    E(a,:)=[];
    W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2+(V(E(:,1),3)-V(E(:,2),3)).^2);
    G=graph(E(:,1),E(:,2),W);
    clear fun
    if length(Source)>10000
        fun=[];
        for s=1:floor(length(Source)/10000)
            dist=distances(G,Source((s-1)*10000+1:s*10000));
            fun=[fun;min(dist,[],1)];
        end
        dist=distances(G,Source(s*10000+1:end));
        fun=[fun;min(dist,[],1)];
        fun=min(fun,[],1);
    else
        dist=distances(G,Source);
        fun=min(dist,[],1);
    end

    fun2=fun*10; %increase the precision
    import edu.stanford.math.plex4.*;
    stream = api.Plex4.createExplicitSimplexStream();
    for i=1:length(V)
        stream.addVertex(i,-fun2(i));
    end
    for i=1:length(E)
        stream.addElement([E(i,1) E(i,2)],max(-fun2(E(i,1)),-fun2(E(i,2))));
    end
    stream.finalizeStream();
    persistence = api.Plex4.getModularSimplicialAlgorithm(1, 2);
    intervals = persistence.computeAnnotatedIntervals(stream);
    gd=intervals.getIntervalsAtDimension(0);
    output_intervals=[file(idx).name(1:end-4),'_',num2str(1),'.txt'];
    diary(output_intervals)
    gd
    diary off
    filename_old=[file(idx).name(1:end-4),'_',num2str(1),'.txt'];
    fileID=fopen(filename_old,'r');
    filename_formatted=[file(idx).name(1:end-4),'_',num2str(1),'_','right_format','.txt'];
    fileID_formatted=fopen(filename_formatted,'w');
    line=fgetl(fileID);
    while ischar(line)
        if length(line)> 1
            line=regexprep(line,'gd','');
            line=regexprep(line,'=', '');
            line=regexprep(line,'\[','');
            line=regexprep(line,'\), ','\n');
            line=regexprep(line,',','');
            line=regexprep(line,'\)\]','');
            line=regexprep(line,'\]','');
            line=regexprep(line,'infinity','0');
            line=regexprep(line,'diary off','');
            if length(line)>0
                fprintf(fileID_formatted,line);
            end
        end
        line=fgetl(fileID);
    end
    fclose(fileID_formatted);
    fclose(fileID);    
end

 %Compute pairwise bottleneck distance BD
 %You can use Javaplex function to compute the bottleneck distance by
 %>> BD(i,j) = edu.stanford.math.plex4.bottleneck.BottleneckDistance.computeBottleneckDistance(gdi,gdj);
 %where gdi is gd for i sample, gdj is gd for j sample
 %But for more efficiently compute it, we computed in Dionysus on our server
 %The instruction can be found https://github.com/danforthcenter/persistent_homology.git

%To compute CC curves
file=dir('diagram*.txt');
 for k=1:length(file)
    diagram=load(file(k).name);
    a=find(diagram(:,1)<-300);
    diagram(a,:)=[];   
    L=diagram(:,2)-diagram(:,1);
    a=find(L>9); %remove any thing less than 1 um
    DG=diagram(a,:);
    y=[];
    for x=-200:10:-10 %level set from 20 um to 1 um
        b=find(DG(:,1)<x&DG(:,2)>x);
        y=[y length(b)];
    end
    BranchFeature=[length(a) sum(L(a))/10 mean(L(a))/10 y];    
    save([file(k).name(9:end-4) '_BranchFeature.mat'],"BranchFeature");
 end 
