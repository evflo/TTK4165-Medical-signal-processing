function s=info2struct(info)
%
% Parse the info string from readecho into a struct
% Torbj?rn Hergum feb. 2010   
% 

    lines=size(info,1);
    mode=0; %0=header, 1=tissue, 2=iq
    s=struct;
	tissue=struct;
	iq=struct;
    for lineNo=1:lines,
        line=info(lineNo,:);
        
        if strfind(line,':'),
            res=regexp(line,':','split');
			res{1}=strtrim(res{1});
			res{2}=strtrim(res{2});
			
                %skip these fields
				if strfind(strtrim(res{1}),'2D tissue parameters'),
					mode=1;
					continue
				elseif strfind(strtrim(res{1}),'IQ parameters'),
					mode=2;
					continue
				end
				if ~isempty(strmatch(strtrim(res{1}),strvcat('Filename','Filelength','Time'))), %#ok<*VCAT>
					continue
				end
				%chopping away whitespace in field name
                res2=regexp(res{1},' ','split');
				if numel(res2)>1,
					res{1}=res2{1};
				end
				%chopping away slashes in field name
                res2=regexp(res{1},'/','split');
				if numel(res2)>1,
					res{1}=res2{1};
				end

					
				
				%Fixing SI units
                res2=regexp(res{2},' ','split');
				if numel(res2)>1,
					if strmatch(res2{2},'MHz'),
						res{2}=num2str(str2num(res2{1})*1e6);
						res{1}=[res{1} '_Hz'];
					end
					if strmatch(res2{2},'fps'),
						res{2}=res2{1};
						res{1}=[res{1} '_fps'];
					end
					if strmatch(res2{2},'m'),
						res{2}=res2{1};
						res{1}=[res{1} '_m'];
					end
					if strmatch(res2{2},'rad'),
						res{2}=res2{1};
						res{1}=[res{1} '_rad'];
					end
					if strmatch(res2{2},'kHz'),
						res{2}=num2str(str2num(res2{1})*1e3);
						res{1}=[res{1} '_Hz'];
					end
					if strmatch(res2{2},'dB'),
						res{2}=res2{1};
						res{1}=[res{1} '_dB'];
					end
					if strmatch(res2{2},'s'),
						res{2}=res2{1};
						res{1}=[res{1} '_s'];
					end
					if strmatch(res2{2},'V'),
						res{2}=res2{1};
						res{1}=[res{1} '_V'];
					end
					if strmatch(res2{2},'='),
						res{2}=res2{3};
					end				
				end
				if strfind(res{1},'StartAngle'),
					%Putting the transducer at the top
					res{2}=num2str(str2num(res{2})-3/2*pi);
				end
				if isempty(str2num(strtrim(res{2}))),
					res{2}=strtrim(res{2});
				else
					res{2}=str2num(strtrim(res{2}));
				end
				
				if mode==0,
					%set names desired names in struct s
					s=setfield(s, res{1}, res{2});
				elseif mode==1,
					tissue=setfield(tissue, strtrim(res{1}), res{2}); %#ok<*SFLD>
					
				elseif mode==2,
					iq=setfield(iq, strtrim(res{1}), res{2});
					
				else
					error('How did you end up here?')
					
				end
        end    
	end
	s.tissue=tissue;
	s.iq=iq;
    
