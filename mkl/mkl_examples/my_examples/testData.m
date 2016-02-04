Points =    [
0.565089 0.975736
0.659576 0.914481
0.064365 0.838268
]; 
gamma = 5;
PointsP = Points./gamma
PointsN = Points./-gamma

Points_e_p = exp(PointsP)
Points_e_n = exp(PointsN)

sum_x_n = 0;
sum_y_n = 0;
sum_x_p = 0;
sum_y_p = 0;


for i=1:1:length(Points_e_n)
    sum_x_n = sum_x_n + Points_e_n(i,1);
    sum_x_p = sum_x_p + Points_e_p(i,1);
    sum_y_n = sum_y_n + Points_e_n(i,2);
    sum_y_p = sum_y_p + Points_e_p(i,2);
end

sum_x_n
sum_y_n
sum_x_p
sum_y_p
    