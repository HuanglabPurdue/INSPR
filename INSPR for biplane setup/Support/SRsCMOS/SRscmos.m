classdef SRscmos < handle
    
    properties (SetAccess = private, GetAccess = public)
        version='7.0';
    end
    
    properties
        %Folder settings
        datafolder;
        resultfolder;
         
        %drift correction
        Cam_pixelsz=120;
        DriftCorrect_frmbin=[];
        DriftCorrect_XYShift=[];
        DriftCorrect_XYZShift=[];
       
    end
    
    properties (SetAccess = public, GetAccess = public)
        frmpfile;
        step_ini=400;
        loc_x;
        loc_y;
        loc_z;
        loc_t;
        loc_bg;
        loc_photons;
        loc_ll;
        loc_crlb;
        loc_x_f;
        loc_y_f;
        loc_z_f;
        loc_t_f;
        loc_step;
        loc_cycle;
        loc_step_f;
        loc_cycle_f;
        st_shifts=[];
        zm = 10;
        Recon_color_highb=3;     
        
    end
   
    
    methods
        function obj=SRscmos(datafolder, resultfolder)
            warning('off')
            %input check
            if nargin<1
                error('Please specify at least: datafolder');
            end
            
            if nargin>2
                error('Please specify: datafolder, resultfolder');
            end
            if nargin == 1
                resultfolder='';
            end
            
            %string check
            if ~ischar(datafolder)&&ischar(resultfolder)
                error('All inputs has to be ''string'' type for: datafolder,resultfolder');
            end
            
            obj.datafolder=datafolder;
            obj.resultfolder=resultfolder;
        end
 
        
        function Perform_DriftCorrection3D(obj)
            x=obj.loc_x.*obj.Cam_pixelsz; %obj.loc_x in pixels
            y=obj.loc_y.*obj.Cam_pixelsz; %obj.loc_y in pixels
            z=obj.loc_z;                  %obj.loc_z in nm
            t=obj.loc_t;
            % xyz is in nm
            reverseflag=0;
            dc_cyc=1;
            unistack=unique(obj.loc_step);
            
            if numel(unistack)>1
                obj.DriftCorrect_frmbin=dc_cyc.*numel(unistack).*obj.frmpfile;
            end
            
            for ss=1:1:numel(unistack)
                display(['Processing optical section number:' num2str(ss)]);
                currst=unistack(ss);
                maskst=(obj.loc_step==currst);
                currx=x(maskst);
                curry=y(maskst);
                currz=z(maskst);
                currt=t(maskst);
                currcy=obj.loc_cycle(maskst);
                [xout,yout,zout,shifts]=driftcorrection_Redun3D(currx,curry,currz,currt,obj.DriftCorrect_frmbin,reverseflag,currcy);
                obj.loc_x_f{ss}=xout;               %obj.loc_x_f in nm
                obj.loc_y_f{ss}=yout;               %obj.loc_y_f in nm
                obj.loc_z_f{ss}=zout;               %obj.loc_z_f in nm
                obj.DriftCorrect_XYZShift{ss}=shifts; % shifts in nm
                
                obj.loc_step_f{ss} = obj.loc_step(maskst); 
                obj.loc_cycle_f{ss} = currcy;
                obj.loc_t_f{ss} = currt;
            end
        end
                
        
        function Perform_stackalignment(obj)
            
            if numel(obj.loc_x_f)==1
                obj.loc_x_f=obj.loc_x_f{1};
                obj.loc_y_f=obj.loc_y_f{1};
                obj.loc_z_f=obj.loc_z_f{1};
            elseif numel(obj.loc_x_f)>1
                pixelsz=50; 
                errorthresh=20;
                cutmeth='nocut';
                iniguess=[0 0 obj.step_ini];
                maskflag=0;
                stksort=[];
                obj.st_shifts=[];
                [shiftx_st,shifty_st,shiftz_st]=stackalignment_Redun3D(stksort,obj.loc_x_f,obj.loc_y_f,obj.loc_z_f,...
                    pixelsz,errorthresh,cutmeth,iniguess,maskflag,[]);

                [xcof]=shiftcoords_stack(obj.loc_x_f,shiftx_st);
                [ycof]=shiftcoords_stack(obj.loc_y_f,shifty_st);
                [zcof]=shiftcoords_stack(obj.loc_z_f,shiftz_st);
                obj.loc_x_f=cat(1,xcof{:});
                obj.loc_y_f=cat(1,ycof{:});
                obj.loc_z_f=cat(1,zcof{:});
                obj.st_shifts=cat(1,shiftx_st,shifty_st,shiftz_st);
            else
                error('No localization detected. Please check detection threshold!');
            end
        end
        
        
    end
end

