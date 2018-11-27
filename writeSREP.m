function writeSREP( filename, nRows, nCols, medialPoints, upPoints, downPoints, crestPoints )
fid = fopen(filename,'wt+');
fprintf(fid, 'model {\n');
fprintf(fid, '\tfigureCount = 1;\n');
fprintf(fid, '\tname = default;\n');
fprintf(fid, '\tfigureTrees {\n');
fprintf(fid, '\t\tcount = 1;\n');
fprintf(fid, '\t\ttree[0] {\n');
fprintf(fid, '\t\t\tattachmentMode = 0;\n');
fprintf(fid, '\t\t\tblendAmount = 0;\n');
fprintf(fid, '\t\t\tblendExtent = 0;\n');
fprintf(fid, '\t\t\tchildCount = 0;\n');
fprintf(fid, '\t\t\tfigureId = 0;\n');
fprintf(fid, '\t\t\tlinkCount = 0;\n');
fprintf(fid, '\t\t}\n');
fprintf(fid, '\t}\n');

fprintf(fid, '\t\tfigure[0] {\n');
fprintf(fid, '\t\tnumRows = %d;\n', nRows);
fprintf(fid, '\t\tnumColumns =  %d;\n', nCols);
fprintf(fid, '\t\tnumLandmarks = 0;\n');
fprintf(fid, '\t\tpositivePolarity = 1;\n');
fprintf(fid, '\t\tpositiveSpace = 1;\n');
fprintf(fid, '\t\tsmoothness = 50;\n');
fprintf(fid, '\t\ttype = QuadFigure;\n');

fprintf(fid, '\t\tcolor {\n');
fprintf(fid, '\t\t\tblue = 0;\n');
fprintf(fid, '\t\t\tgreen = 0.5;\n');
fprintf(fid, '\t\t\tred = 0;\n');
fprintf(fid, '\t\t}\n');

for i = 1:nRows
    for j = 1:nCols
        fprintf(fid,'\t\tprimitive[%d][%d]{\n',i-1,j-1);
        fprintf(fid,'\t\t\tselected = 1;\n');
        if i == 1 || i == nRows || j == 1 || j == nCols
            fprintf(fid, '\t\t\ttype = EndPrimitive;\n');
        else
            fprintf(fid, '\t\t\ttype = StandardPrimitive;\n');
        end
        current_medial_x = medialPoints(i,j,1);
        current_medial_y = medialPoints(i,j,2);
        current_medial_z = medialPoints(i,j,3);
        medial_v = [current_medial_x, current_medial_y, current_medial_z];
        fprintf(fid,'\t\t\tx = %f;\n',current_medial_x);
        fprintf(fid,'\t\t\ty = %f;\n',current_medial_y);
        fprintf(fid,'\t\t\tz = %f;\n',current_medial_z);

        current_up_x = upPoints(i,j,1);
        current_up_y = upPoints(i,j,2);
        current_up_z = upPoints(i,j,3);
        up_v = [current_up_x, current_up_y, current_up_z];
        up_direction_v = normr(up_v - medial_v);
        up_length = norm(up_v - medial_v);
        
        
        current_down_x = downPoints(i,j,1);
        current_down_y = downPoints(i,j,2);
        current_down_z = downPoints(i,j,3);
        down_v = [current_down_x, current_down_y, current_down_z];
        down_direction_v = normr(down_v - medial_v);
        down_length = norm(down_v - medial_v);
        
        current_crest_x = crestPoints(i,j,1);
        current_crest_y = crestPoints(i,j,2);
        current_crest_z = crestPoints(i,j,3);
        crest_v = [current_crest_x, current_crest_y, current_crest_z];
        crest_direction_v = normr(crest_v - medial_v);
        crest_length = norm(crest_v - medial_v);
        if i == 1 || i == nRows || j == 1 || j == nCols
            if crest_length < up_length || crest_length < down_length
                warning(sprintf('At (%d, %d), the crest length is not the max\n', i,j));
                crest_length = max([crest_length, up_length, down_length]);
            end
        end
        fprintf(fid, '\t\t\tux[0] = %f;\n', up_direction_v(1));
        fprintf(fid, '\t\t\tux[1] = %f;\n', down_direction_v(1));
        fprintf(fid, '\t\t\tux[2] = %f;\n', crest_direction_v(1));        
        
        fprintf(fid, '\t\t\tuy[0] = %f;\n', up_direction_v(2));
        fprintf(fid, '\t\t\tuy[1] = %f;\n', down_direction_v(2));
        fprintf(fid, '\t\t\tuy[2] = %f;\n', crest_direction_v(2));
                
        fprintf(fid, '\t\t\tuz[0] = %f;\n', up_direction_v(3));
        fprintf(fid, '\t\t\tuz[1] = %f;\n', down_direction_v(3));
        fprintf(fid, '\t\t\tuz[2] = %f;\n', crest_direction_v(3));
        
        fprintf(fid, '\t\t\tr[0] = %f;\n', up_length);
        fprintf(fid, '\t\t\tr[1] = %f;\n', down_length);
%         
        if i == 1 || i == nRows || j == 1 || j == nCols
            fprintf(fid, '\t\t\tr[2] = %f;\n', crest_length);
        end
        fprintf(fid, '\t\t}\n');
    end
end
fprintf(fid, '\t}\n');
fprintf(fid, '\ttransformation {\n');
fprintf(fid, '\t\tscale = 1;\n');
fprintf(fid, '\t\trotation {\n');
fprintf(fid,'\t\t\tw = 1;\n');
fprintf(fid,'\t\t\tx = 0;\n');
fprintf(fid,'\t\t\ty = 0;\n');
fprintf(fid,'\t\t\tz = 0;\n');
fprintf(fid,'\t\t}\n');
fprintf(fid,'\t\ttranslation {\n');
fprintf(fid,'\t\t\tx = 0;\n');
fprintf(fid,'\t\t\ty = 0;\n');
fprintf(fid,'\t\t\tz = 0;\n');
fprintf(fid,'\t\t}\n');
fprintf(fid,'\t}\n');
fprintf(fid,'}\n');
fclose(fid);
end

