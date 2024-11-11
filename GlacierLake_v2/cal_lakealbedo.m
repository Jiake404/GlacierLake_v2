function [lake_albedo] = cal_lakealbedo(albedo_num,lake_depth)
if albedo_num == 0
    lake_albedo = (9702 + 1000.*exp(3.6*lake_depth))...
            /(-539 + 20000.*exp(3.6*lake_depth)); 
else
lake_albedo = 0.1820+0.3512.*exp(-1.0774.*lake_depth);
end
