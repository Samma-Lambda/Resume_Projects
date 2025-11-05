#constraint 36 doesn't incorporate drone weights into the weight capacity of truck


function [solution_vector,customers, truck_set,drone_set, V_L, V_R, C, x, s_U, r, y, t_TA,t_TL,t_UA,t_UL,w_i,Q,q,e,s_T,w_U] = LP_truck_drone()
  rand("seed", 123)

  Verbose = true;
  trucks = 2;
  drones = 2;
  customers = 10;
  truck_speeds = [50,50];
  drone_speeds = [80,70];
  BRT = 2;
  M = 10000;
  epsilon = 0.00001;
  V = 1:(customers + 2);
  V_L = V(1:end-1);
  V_R = V(2:end);
  drone_set = 1:drones;
  truck_set = 1:trucks;
  Z = 1200 * ones(1:trucks);
  Q = [3,4];
  e = [0.623,0.844];

  #Generate customer locations
  width = 5;
  height = 5;
  [customer_points,distance_matrix] = make_distance_matrix(customers,width,height);


  #generate customer types
  [C,C_T, C_U, C_F] = customer_sets(customers);


  #generate weights
  w_i = Generate_Weight(C,C_T,C_U,C_F,Q,Z,trucks,drones);

  if Verbose == true
    disp("customer points");
    disp(customer_points);
    disp("customer weights");
    disp(w_i);
    disp("All Customers");
    disp(C);
    disp("Truck Only Customers:");
    disp(C_T);
    disp("Drone Only Customers:");
    disp(C_U);
    disp("Flexible Customers:");
    disp(C_F);
    disp("Customer Weights")
    disp(w_i)
    disp("V: ")
    disp(V)
    disp("V_L:")
    disp(V_L)
    disp("V_R")
    disp(V_R)
    disp("Drone set")
    disp(drone_set)
    disp("Drone Speeds")
    disp(drone_speeds)
    disp("Truck set ")
    disp(truck_set)
    disp("Truck speeds");
    disp(truck_speeds)
    disp("Battery Replacement Time: ")
    disp(BRT)
    disp("truck capacities:")
    disp(Z)
    disp("Drone Weight Capacities")
    disp(Q)
    disp("Drone energy capacities")
    disp(e)
  endif



  t_T = [];
  for i=V_L
    for j=V_R
       for t = truck_set
          t_T(i,j,t) = (distance_matrix(i,j)/truck_speeds(t))*60;
       endfor
    endfor
  endfor

  t_U = [];
  for i=V_L
    for j=V_R
       for o = drone_set
          t_U(i,j,o) = (distance_matrix(i,j)/drone_speeds(o))*60;
       endfor
    endfor
  endfor

  sigma_t = [];
  for t = truck_set
    for i = C
       sigma_t(i,t)=5;
    endfor
  endfor

  sigma_u = [];
  for o = drone_set
    for i = C
       sigma_u(i,o)=3;
    endfor
  endfor


  #######################
  #---------DV's--------#
  #######################
  x = [];
  counter = 1;
  variable_types = char();

  for t = truck_set
    for i = V_L
      for j = V_R
        if i!=j
          x(i,j,t)=counter;
          variable_types = [variable_types; "B"];
          counter = counter + 1;
        endif
      endfor
    endfor
  endfor
  disp("x");
  disp(counter);

  y = [];
  for o = drone_set
    for t = truck_set
      for i  = V_L
        for j = V_R
          if i != j
            y(i,j,o,t)= counter;
            variable_types = [variable_types; "B"];
            counter = counter + 1;
          endif
        endfor
      endfor
    endfor
  endfor
  disp("y");
  disp(counter);

  s_T = [];
  for t = truck_set
    for i = C
      s_T(i,t)= counter;
      variable_types = [variable_types; "B"];
      counter = counter + 1;
    endfor
  endfor
  disp("S_T");
  disp(counter);

  s_U = [];
  for o = drone_set
    for t = truck_set
      for j = V_L
        for i = C
          if i != j
            s_U(i,j,o,t) = counter;
            variable_types = [variable_types; "B"];
            counter =  counter + 1;
          endif
        endfor
      endfor
    endfor
  endfor
  disp("S_U");
  disp(size(s_U));
  disp(counter);

  r = [];
  for o = drone_set
    for t = truck_set
      for i = V_L
        for j = V_R
          if i != j
             r(i,j,o,t) = counter;
             variable_types = [variable_types; "B"];
             counter = counter + 1;
          endif
        endfor
      endfor
    endfor
  endfor
  disp("r")
  disp(counter)

  t_UA = [];
  for o = drone_set
    for t = truck_set
      for i = V_R
        t_UA(i,o,t)=counter;
        variable_types = [variable_types; "C"];
        counter = counter +1;
      endfor
    endfor
  endfor
  disp("t_UA")
  disp(counter)


  t_UL = [];
  for o = drone_set
    for i = V_L
      for t = truck_set
        t_UL(i,o,t)=counter;
        variable_types = [variable_types; "C"];
        counter = counter + 1;
      endfor
    endfor
  endfor
  disp("t_UL")
  disp(counter)

  t_TA = [];
  for t = truck_set
    for i = V_R
      t_TA(i,t) = counter;
      variable_types = [variable_types; "C"];
      counter = counter + 1;
    endfor
  endfor
  disp("t_TA")
  disp(counter)

  t_TL = [];
  for t = truck_set
    for i = V_L
      t_TL(i,t) = counter;
      variable_types = [variable_types; "C"];
      counter = counter + 1;
    endfor
  endfor
  disp("t_TL")
  disp(counter)

  w_U = [];
  for i = V
    for o = drone_set
      for t = truck_set
        w_U(i,o,t)=counter;
        variable_types = [variable_types; "C"];
        counter = counter + 1;
      endfor
    endfor
  endfor
  disp("w_U")
  disp(counter)

  q = [];
  for o = drone_set
    for t = truck_set
      for i = V_L
        q(t,i,o)=counter;
        variable_types = [variable_types; "C"];
        counter = counter + 1;
      endfor
    endfor
  endfor
  disp("q")
  disp(counter)

  u = [];
  for t = truck_set
    for i = V
      u(t,i)= counter;
      variable_types = [variable_types; "I"];
      counter = counter + 1;
    endfor
  endfor
  disp("u")
  disp(counter)

  p_t = [];
  for t =  truck_set
    for i = V_L
      for j = V_R
        if i != j
           p_t(t,i,j)=counter;
           variable_types = [variable_types; "B"];
           counter = counter + 1;
        endif
      endfor
    endfor
  endfor
  disp("p")
  disp(counter)

  big_T = counter;
  variable_types = [variable_types; "C"];
  disp(counter);


  ###############################
  #---------Constraints---------#
  ###############################

  constraint_matrix = [];
  RHS = [];
  Constraint_Types = char();

  #our objective function
  c = zeros(1, big_T);
  c(big_T) = 1;

  #minimize drone weights
  for i = V_L
    for o = drone_set
      for t = truck_set
        c(w_U(i,o,t))=epsilon;
      endfor
    endfor
  endfor

  #minimize truck arrival/depature times
  for t = truck_set
    for i = V_R
      c(t_TA(i,t))=epsilon;
    endfor
  endfor

  for t = truck_set
    for i = V_L
      c(t_TL(i,t)) = epsilon;
    endfor
  endfor

  #minimize drone arrival/depature times
  for o = drone_set
    for t = truck_set
      for i = V_R
        c(t_UA(i,o,t))=epsilon;
      endfor
    endfor
  endfor

  for o = drone_set
    for i = V_L
      for t = truck_set
        c(t_UL(i,o,t))=epsilon;
      endfor
    endfor
  endfor



S_counter = 0;
U_counter = 0;

#Constraint 1
  for t = truck_set
      row = zeros(1,counter);
      row(t_TA(customers+2,t))=1;
      row(big_T)=-1;

      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 0];
      Constraint_Types = [Constraint_Types; 'U'];
      U_counter = U_counter + 1;
  endfor
disp("1");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);

#Constraint 2
  for t = truck_set
    for o = drone_set
      row = zeros(1,counter);
      row(t_UA(customers+2,o,t))=1;
      row(big_T)=-1;
      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 0];
      Constraint_Types = [Constraint_Types; 'U'];
      U_counter = U_counter + 1;
    endfor
  endfor
disp("2");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);


%%%%%%  #constraint 3
num_C = length(C);
constraint_matrix_temp = zeros(num_C, counter);
RHS_temp = ones(num_C, 1); % All RHS values are 1
Constraint_Types_temp = repmat('S', num_C, 1); % All 'S'

for idx = 1:num_C
    i = C(idx);
    row = zeros(1, counter);

    for o = drone_set
        for t = truck_set
            for j = V_L
                if i ~= j
                    row(s_U(i, j, o, t)) = 1;
                end
            end
        end
    end

    for t = truck_set
        row(s_T(i, t)) = 1;
    end

    S_counter = S_counter + 1;
    constraint_matrix_temp(idx, :) = row;
end
%%
%%% Append preallocated arrays
constraint_matrix = [constraint_matrix; constraint_matrix_temp];
RHS = [RHS; RHS_temp];
Constraint_Types = [Constraint_Types; Constraint_Types_temp];

disp("3");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);
%%%%
%%%%%%
%%%%%%%%%%
%%%%%%  #constraint 4
  for i = C_T
     row = zeros(1, counter);
     for t = truck_set
        row(s_T(i,t))=1;
     endfor
     constraint_matrix = [constraint_matrix; row];
     RHS = [RHS; 1];
     Constraint_Types = [Constraint_Types; 'S'];
     S_counter = S_counter + 1;
  endfor
  disp("4");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%
%%%%  #constraint 5
  for i = C_U
    row = zeros(1, counter);
    for o = drone_set
      for t = truck_set
        for j = V_L
          if i!=j
             row(s_U(i,j,o,t))=1;
          endif
        endfor
      endfor
    endfor
    constraint_matrix = [constraint_matrix; row];
    RHS = [RHS; 1];
    Constraint_Types = [Constraint_Types; 'S'];
    S_counter = S_counter + 1;
  endfor
  disp("5");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%
%%%%%%%%  #constraint 6
  for t = truck_set
      row = zeros(1, counter);
      for i = V_R
          row(x(1,i,t))=1;
      endfor
      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 1];
      Constraint_Types = [Constraint_Types; 'S'];
      S_counter = S_counter + 1;
  endfor
  disp("6");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%%%  #constraint 7
  for t = truck_set
     row = zeros(1,counter);
     for i = V_L
       row(x(i,customers+2,t))=1;
     endfor
     constraint_matrix = [constraint_matrix; row];
     RHS = [RHS; 1];
     Constraint_Types = [Constraint_Types; 'S'];
     S_counter = S_counter + 1;
  endfor

  disp("7");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%
%%%%%%  #constraint 8 A
  for j = C
    for t = truck_set
      row = zeros(1,counter);
      for i = V_L
        if i!=j
          row(x(i,j,t))=1;
        endif
      endfor

      for h = V_R
        if h!=j
          row(x(j,h,t))=-1;
        endif
      endfor
      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 0];
      Constraint_Types = [Constraint_Types; 'S'];
      S_counter = S_counter + 1;
    endfor
  endfor

  disp("8 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%  #constraint 8 B
  for j = C
     for t = truck_set
        row=zeros(1,counter);
        for i = V_L
           if i != j
             row(x(i,j,t))=1;
           endif
        endfor
        constraint_matrix = [constraint_matrix; row];
        RHS = [RHS; 1];
        Constraint_Types = [Constraint_Types; 'U'];
        U_counter = U_counter + 1;
     endfor
  endfor

  disp("8 B");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%  #constraint 9
  for i = V_L
    for j = V_R
      if i != j
         for t = truck_set
           row = zeros(1, counter);
           row(u(t,i)) = 1;
           row(u(t,j))=-1;
           row(x(i,j,t))= M;
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; M - 1];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
         endfor
      endif
    endfor
  endfor
  disp("9");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%
%%%%%%  #constraint 10
  for i = V_L
    for j = V_R
      if i != j
        for t = truck_set
          row = zeros(1,counter);
          row(u(t,i))=-1;
          row(u(t,j))=1;
          row(p_t(t,i,j))= -M;
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; -1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endif
    endfor
  endfor
  disp("10");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%%%%%  #constraint 11
  for i = V_L
    for j = V_R
      if i != j
        for t = truck_set
          row = zeros(1, counter);
          row(u(t,i)) = 1;
          row(u(t,j)) = -1;
          row(p_t(t,i,j)) = M;
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; M-1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endif
    endfor
  endfor
  disp("11");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%%%%%
%%%%#constraint 12 A
  for t = truck_set
    for i = C
      for j = C
        if i!=j
           row = zeros(1,counter);
           row(p_t(t,i,j))=1;
           row(p_t(t,j,i))=1;
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 1];
           Constraint_Types = [Constraint_Types; 'S'];
           S_counter = S_counter + 1;
        endif
      endfor
    endfor
  endfor

  disp("12 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%  #constraint 12 B
  for t = truck_set
      for j = V_R
        if i!=j
           row = zeros(1,counter);
           row(p_t(t,1,j))=1;
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 1];
           Constraint_Types = [Constraint_Types; 'S'];
           S_counter = S_counter + 1;
        endif
      endfor
  endfor
  disp("12 B");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);

  #constraint 12 C
  for t = truck_set
    for i = V_L
        if i!=j
           row = zeros(1,counter);
           row(p_t(t,i,customers+2))=1;
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 1];
           Constraint_Types = [Constraint_Types; 'S'];
           S_counter = S_counter + 1;
        endif
    endfor
  endfor
  disp("12 C");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%%%  #constraint 13
  for o = drone_set
    for t = truck_set
      for i = C
        row = zeros(1,counter);
        for p  = V_L
          if p != i
            row(y(p,i,o,t)) = -1;
          endif
        endfor

        for j = V_L
           if i != j
              row(s_U(i,j,o,t))=1;
           endif
        endfor
        constraint_matrix = [constraint_matrix; row];
        RHS = [RHS; 0];
        Constraint_Types = [Constraint_Types; 'U'];
        U_counter = U_counter + 1;
      endfor
    endfor
  endfor
  disp("13");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%
%%%%%%%%
%%%%%%
%%%%  #constraint 14
%%%%%%
  for i = C
    for t = truck_set
      row =  zeros(1,counter);
      row(s_T(i,t))=-1;
      for j = V_L
        if i != j
          row(x(j,i,t))=1;
        endif
      endfor
      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 0];
      Constraint_Types = [Constraint_Types; 'S'];
      S_counter = S_counter + 1;
    endfor
  endfor
  disp("14");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%
  #constraint 15
  for i = V_L
     for o = drone_set
        for t = truck_set
           row = zeros(1,counter);
           for j = V_R
              if i != j
                 row(r(i,j,o,t))=1;
              endif
           endfor
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 1];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
        endfor
     endfor
   endfor
  disp("15");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%   #constraint 16
   for i = V_R
     for o = drone_set
        for t = truck_set
           row = zeros(1,counter);
           for j = V_L
              if i != j
                 row(r(j,i,o,t))=1;
              endif
           endfor
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 1];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
        endfor
     endfor
   endfor
  disp("16");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%## Constraint 17 A
for p = C
  for o = drone_set
     for t = truck_set
       row = zeros(1,counter);
       for i = V_L
         if i!=p
            row(y(i,p,o,t))=-1;
         endif
       endfor

       for j = V_R
          if p!=j
            row(y(p,j,o,t))=1;
          endif
       endfor

       row(s_T(p,t))=-M;
       constraint_matrix = [constraint_matrix; row];
       RHS = [RHS; 0];
       Constraint_Types = [Constraint_Types; 'U'];
       U_counter = U_counter + 1;
     endfor
  endfor
endfor
disp("17 A");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);

%%%%%%
## Constraint 17 B NEW
for p = C
  for o = drone_set
     for t = truck_set
       row = zeros(1,counter);
       for i = V_L
         if i!=p
            row(y(i,p,o,t))=1;
         endif
       endfor


       for j = V_R
         if p!=j
            row(y(p,j,o,t))=-1;
         endif
       endfor

       row(s_T(p,t))=-M;


       constraint_matrix = [constraint_matrix; row];
       RHS = [RHS; 1];
       Constraint_Types = [Constraint_Types; 'U'];
       U_counter = U_counter + 1;
     endfor
  endfor
endfor
disp("17 B");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);

%%
%%%%## Constraint 17 *
for p = C
  for o = drone_set
     for t = truck_set
       row = zeros(1,counter);
       for i = V_L
          if i!=p
            row(y(i,p,o,t))=-1;
          endif
       endfor

       for j = V_R
          if j!=p
             row(y(p,j,o,t))=1;
          endif
       endfor

       row(s_T(p,t))=-M;

       constraint_matrix = [constraint_matrix; row];
       RHS = [RHS; 0];
       Constraint_Types = [Constraint_Types; 'U'];
       U_counter = U_counter + 1;
     endfor
  endfor
endfor
disp("17 *");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%## Constraint 17 **
for p = C
  for o = drone_set
     for t = truck_set
       row = zeros(1,counter);
       for i = V_L
          if i!=p
            row(y(i,p,o,t))=1;
          endif
       endfor

       for j = V_R
          if j!=p
             row(y(p,j,o,t))=-1;
          endif
       endfor

       row(s_T(p,t))=-M;

       constraint_matrix = [constraint_matrix; row];
       RHS = [RHS; 0];
       Constraint_Types = [Constraint_Types; 'U'];
       U_counter = U_counter + 1;
     endfor
  endfor
endfor
disp("17 **");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);
%%%%
%%%%
%%%%%%
%%%%%%%%   #constraint 18
   for i = C
     for j = C
       if i != j
       row = zeros(1,counter);
       for o = drone_set
         for t = truck_set
            if i!=j
              row(y(i,j,o,t))=1;
            endif
         endfor
       endfor
       for o = drone_set
         for t = truck_set
            if i != j
              row(y(j,i,o,t))=1;
            endif
         endfor
       endfor
       constraint_matrix = [constraint_matrix; row];
       RHS = [RHS; 1];
       Constraint_Types = [Constraint_Types; 'U'];
       U_counter = U_counter + 1;
       endif
     endfor
   endfor
  disp("18");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%%%   #constraint 19 A
   for i = C
      for j = C
        if i != j
        for o = drone_set
           for t = truck_set
             row = zeros(1,counter);

             if j != i
               row(s_U(j,i,o,t))=-1;
             endif

             if j != i
               row(y(i,j,o,t))=M;
             endif

             row(s_T(i,t))=M;

             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; (2*M)-1];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
           endfor
        endfor
        endif
      endfor
   endfor
  disp("19 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%
%%  #constraint 19 B
  for j = C
     for o = drone_set
       for t = truck_set
          row = zeros(1,counter);
          row(s_U(j,1,o,t))=-1;
          row(y(1,j,o,t))=M;
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; M-1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
         endfor
      endfor
   endfor
  disp("19 B");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%   #constraint 20 A
   for i = C
    for o = drone_set
      for t = truck_set
         row= zeros(1,counter);
         for l = V_R
           if i!=l
             row(r(i,l,o,t))=-1;
           endif
         endfor

         for j = C
           if i!=j
             row(y(i,j,o,t))=M;
           endif
         endfor
         row(s_T(i,t))= M;

         constraint_matrix = [constraint_matrix; row];
         RHS = [RHS; (2*M)-1];
         Constraint_Types = [Constraint_Types; 'U'];
         U_counter = U_counter + 1;
      endfor
    endfor
   endfor
  disp("20 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);

  #constraint 20 B
    for o = drone_set
      for t = truck_set
         row= zeros(1,counter);
         for l = V_R
             row(r(1,l,o,t))=-1;
         endfor

         for j = C
           row(y(1,j,o,t))=M;
         endfor

         constraint_matrix = [constraint_matrix; row];
         RHS = [RHS; M-1];
         Constraint_Types = [Constraint_Types; 'U'];
         U_counter = U_counter + 1;
      endfor
   endfor
   disp("20 B");
   printf("Number of 'S': %d\n", S_counter);
   printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%%%%%    #constraint 21 A
    for j = C
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          for l = V_L
            if l !=j
              row(r(l,j,o,t))=-1;
            endif
          endfor

          for i = C
            if i!=j
              row(y(i,j,o,t))=M;
            endif
          endfor

          row(s_T(j,t))=M;

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; (2*M)-1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
  disp("21 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%
    #21 B
    for o = drone_set
      for t = truck_set
        row = zeros(1,counter);
        for l = V_L
          row(r(l,customers+2,o,t))=-1;
        endfor

        for i = C
          row(y(i,customers+2,o,t))=M;
        endfor


        constraint_matrix = [constraint_matrix; row];
        RHS = [RHS; M-1];
        Constraint_Types = [Constraint_Types; 'U'];
        U_counter = U_counter + 1;
      endfor
   endfor
  disp("21 B");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);

  #constraint 22 A
    for i = V_L
      for j = C
      if i!=j
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          for j = C
            if i !=j
              row(s_U(j,i,o,t))=1;
            endif
          endfor

          for l = V_R
            if i!=l
              row(r(i,l,o,t))=-M;
            endif
          endfor

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; 0];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
      endif
      endfor
    endfor
  disp("22 A");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%  #constraint 22B
  for o = drone_set
    for t = truck_set
       for i = V_L
         for j = C
            if i!=j
               row = zeros(1,counter);
               row(t_UL(i,o,t))=1;
               row(t_UA(j,o,t))=-1;
               row(s_U(j,i,o,t))=M;

               constraint_matrix = [constraint_matrix; row];
               RHS = [RHS; M];
               Constraint_Types = [Constraint_Types; 'U'];
               U_counter = U_counter + 1;
            endif
         endfor
       endfor
    endfor
  endfor
  disp("22 B");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);
%%%%%%%%
%%%%%%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    #constraint 23  ##Check
    for i = V_L
       for o = drone_set
          for t = truck_set
             row = zeros(1,counter);
             for j = C
                if i!=j
                  row(s_U(j,i,o,t))= -1;
                endif
             endfor

             for l = V_R
                if i != l
                  row(r(i,l,o,t))=1;
                endif
             endfor

             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; 0];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
          endfor
       endfor
    endfor
    disp("23");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%
%%%%%%%%%%%%    #constraint 24  ##Check
    for i = V_L
      for j =V_R
        if i != j
        for t = truck_set
           for o = drone_set
             row = zeros(1,counter);
             row(u(t,j))=-1;
             row(u(t,i))=1;
             if i != j
               row(r(i,j,o,t))=M;
             endif
             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; M-1];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
           endfor
        endfor
        endif
      endfor
    endfor
    disp("24");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%    #constraint 25  ##Check
    for i = V_L
       for t = truck_set
         for o = drone_set
           row = zeros(1,counter);
           for j = V_R
              if i != j
                row(x(i,j,t))=-1;
              endif
           endfor
           for l = V_R
             if l!=i
               row(r(i,l,o,t))= M;
             endif
           endfor
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; M-1];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
         endfor
       endfor
    endfor
    disp("25");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%
%%%%%%%%%%%%%%%%    #constraint 26   ##Check
    for i = V_R
      for t = truck_set
         for o = drone_set
            row = zeros(1,counter);
            for j = V_L
              if i !=j
                row(x(j,i,t))=-1;
              endif
            endfor

            for l = V_L
              if l != i
                row(r(l,i,o,t))=M;
              endif
            endfor
            constraint_matrix = [constraint_matrix; row];
            RHS = [RHS; M-1];
            Constraint_Types = [Constraint_Types; 'U'];
            U_counter = U_counter + 1;
         endfor
      endfor
    endfor
    disp("26");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
% Constraint 27: Adds time to trucks when replacing batteries
for i = V_L
  for p = C
    if i != p
      for j = V_R
        if i != j
          for o = drone_set
            for t = truck_set

              row = zeros(1, counter);

              % launch/retrieve terms
              row(t_UL(p, o, t)) = -1;
              row(t_UA(j, o, t)) =  1;

              % battery-return-time terms
              for l = drone_set
                for m = C
                  if m != j
                    row(y(m, j, l, t)) = BRT;
                  end
                end
              end

              % big-M terms
              row(r(i, j, o, t)) = M;

              for q_t = V_R
                if p != q_t
                  row(r(p, q_t, o, t)) = row(r(p, q_t, o, t)) + M;
                end
              end

              row(p_t(t, i, p)) = M;

              % Append to global constraint structures
              constraint_matrix = [constraint_matrix; row];
              RHS               = [RHS; 3 * M];
              Constraint_Types  = [Constraint_Types; 'U'];
              U_counter         = U_counter + 1;

            end
          end
        end
      end
    end
  end
end

disp("27");
printf("Number of 'S': %d\n", S_counter);
printf("Number of 'U': %d\n", U_counter);





%%%%%%%%%%%%%%
%%%%%%%%%%    #constraint 28   ##Make truck arrival time work
    for i = V_L
      for j = V_R
         if i!=j
         for t = truck_set
         row = zeros(1,counter);
         row(t_TA(j,t))=-1;
         row(t_TL(i,t))=1;

         row(x(i,j,t))=M;


         constraint_matrix = [constraint_matrix; row];
         RHS = [RHS; M - t_T(i,j,t)];
         Constraint_Types = [Constraint_Types; 'U'];
         U_counter = U_counter + 1;
         endfor
         endif
      endfor
    endfor
    disp("28");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

    #constraint 29 A
    for i = C
      for t = truck_set
        row=zeros(1,counter);
        row(t_TL(i,t))=-1;
        row(t_TA(i,t))=1;

        row(s_T(i,t))=M;

        for o = drone_set
           for j = C
             if i!=j
               row(y(j,i,o,t))=BRT;
             endif
           endfor
        endfor

        constraint_matrix = [constraint_matrix; row];
        RHS = [RHS; M - sigma_t(i,t)];
        Constraint_Types = [Constraint_Types; 'U'];
        U_counter = U_counter + 1;
      endfor
    endfor
    disp("29 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

    #constraint 29 B
    for i = C
      for o = drone_set
         for t = truck_set
            row = zeros(1,counter);

            row(t_TL(i,t))=-1;
            row(t_UA(i,o,t))=1;

            for j = C
              if i!=j
                row(y(j,i,o,t))=BRT;
              endif
            endfor

            row(s_T(i,t))=M;

            constraint_matrix = [constraint_matrix; row];
            RHS = [RHS; M ];
            Constraint_Types = [Constraint_Types; 'U'];
            U_counter = U_counter + 1;

         endfor
      endfor
    endfor
    disp("29 B");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);


#constraint 30  ##Check
    for i = V_L
      for j = V_R
        if i!=j
        for o = drone_set
          for t = truck_set
          row = zeros(1,counter);
          row(t_UA(j,o,t))=-1;
          row(t_UL(i,o,t))=1;
          row(y(i,j,o,t))=M;

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; M - t_U(i,j,o)];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
          endfor
        endfor
        endif
      endfor
    endfor
    disp("30");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%    #constraint 31  ##Check
    for i = C
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          row(t_UL(i,o,t))=-1;
          row(t_UA(i,o,t))=1;

          for j = V_L
            if i != j
              row(s_U(i,j,o,t))=M;
            endif
          endfor
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; M - sigma_u(i,o)];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("31 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

#constraint 31 B
     for i = C
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          row(t_UL(i,o,t))=-1;
          row(t_UA(i,o,t))=1;

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; 0];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("31 B");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%%%    #constraint 32 ##Check
    for i = V_L
      for t = truck_set
        for o = drone_set
          row = zeros(1,counter);
          row(t_UL(i,o,t))=-1;
          row(t_TL(i,t))=1;
          for l = V_R
            if i != l
              row(r(i,l,o,t))=M;
            endif
          endfor
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; M];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("32");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%
%%%%    #constraint 33 ##Check
    for i = V_L
      for t = truck_set
        for o = drone_set
           row =  zeros(1,counter);
           row(t_UL(i,o,t))=1;
           row(t_TL(i,t))=-1;
           for l = V_R
             if i!=l
               row(r(i,l,o,t))=M;
             endif
           endfor
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; M];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("33");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%    #constraint 34 ##CHECK
    for i = C
       for t = truck_set
         for o = drone_set
           row = zeros(1,counter);
           row(t_TL(i,t))=-1;
           row(t_UA(i,o,t))=1;
           row(s_T(i,t))=M;
           for j = C
             if j != i
               row(y(j,i,o,t))=M;
             endif
           endfor
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; M*2];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
         endfor
       endfor
    endfor
    disp("34");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);


%%%%%%
%%%%%%    #constraint 35   ##truck returns back after 450 minutes
    for t = truck_set
      row = zeros(1,counter);
      row(t_TA(customers+2,t))=1;

      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; 450];
      Constraint_Types = [Constraint_Types; 'U'];
      U_counter = U_counter + 1;
    endfor
    disp("35");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%
%%%%    #constraint 36 weight of all drone launches is less than drone capacity
    for j = V_L
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          for i = C
            if i != j
              row(s_U(i,j,o,t))=w_i(i);
            endif
          endfor
          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; Q(o)];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("36");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

 #constraint 37    the weight of all drone launches and truck deliveries is less than the capacity of the truck
    for t = truck_set
      row = zeros(1,counter);
      for o = drone_set
        for i = C
           for j = V_L
             if i != j
                row(s_U(i,j,o,t))=w_i(j);
             endif
           endfor
        endfor
      endfor
      for i = C
        row(s_T(i,t))=w_i(i);
      endfor
      constraint_matrix = [constraint_matrix; row];
      RHS = [RHS; Z(t)];
      Constraint_Types = [Constraint_Types; 'U'];
      U_counter = U_counter + 1;
    endfor
    disp("37");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
    #constraint 38 A
    for i = V_L
      for j = C
        if i!=j
        for o = drone_set
           for t = truck_set
             row = zeros(1,counter);
             row(w_U(j,o,t))=-1;

             for m = C
               if m!=i
                  row(s_U(m,i,o,t))= w_i(m);
               endif
             endfor

             if i!=j
              row(y(i,j,o,t))=M;
             endif

             for l = V_R
              if i!=l
                row(r(i,l,o,t))=M;
              endif
             endfor

             if j != i
                row(s_U(j,i,o,t))=row(s_U(j,i,o,t))+M;
             endif

             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; (3*M)+w_i(j)];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
           endfor
        endfor
        endif
      endfor
    endfor
    disp("38 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

%%  # constraint 38 C
  for i = V_L
     for o = drone_set
       for t = truck_set
         row = zeros(1,counter);
         row(w_U(i,o,t))=-1;
         for m = C
            if i != m
              row(s_U(m,i,o,t))=w_i(m);
            endif
         endfor

         constraint_matrix = [constraint_matrix; row];
         RHS = [RHS; 0];
         Constraint_Types = [Constraint_Types; 'U'];
         U_counter = U_counter + 1;
       endfor
     endfor
  endfor
  disp("38 C");
  printf("Number of 'S': %d\n", S_counter);
  printf("Number of 'U': %d\n", U_counter);


#constraint 39 A
    for i = C
      for j = C
         if i!=j
         for l = V_L
            for o = drone_set
              for t = truck_set
              row = zeros(1, counter);
              row(w_U(j,o,t))=-1;
              row(w_U(i,o,t))=1;
              if i != j
                row(y(i,j,o,t))=M;
              endif
              if j != l
                 row(s_U(j,l,o,t))=row(s_U(j,l,o,t)) + M;
              endif
              if i !=l
                row(s_U(i,l,o,t))=row(s_U(i,l,o,t)) + M;
              endif

              constraint_matrix = [constraint_matrix; row];
              RHS = [RHS; (3*M)+w_i(j)];
              Constraint_Types = [Constraint_Types; 'U'];
              U_counter = U_counter + 1;
              endfor
            endfor
         endfor
         endif
      endfor
    endfor
    disp("39 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);






%%%% #constraint 40 A #check
    for i = C
       for o = drone_set
          for t = truck_set
             row = zeros(1,counter);
             row(q(t,i,o)) = -1;
             row(s_T(i,t))=M;
             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; M-e(o)];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
          endfor
       endfor
    endfor
    disp("40 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%
%%# constraint 40 B #Check
 for o = drone_set
   for t = truck_set
     row = zeros(1,counter);
     row(q(t,1,o))=-1;

     constraint_matrix = [constraint_matrix; row];
     RHS = [RHS; -e(o)];
     Constraint_Types = [Constraint_Types; 'U'];
     U_counter = U_counter + 1;
   endfor
 endfor
 disp("40 B");
 printf("Number of 'S': %d\n", S_counter);
 printf("Number of 'U': %d\n", U_counter);
#constraint 41 A #check
    for i = V_L
      for j = C
        if i!=j
        for o = drone_set
          for t = truck_set
            row = zeros(1,counter);
            row(q(t,j,o))=-1;
            row(q(t,i,o))=1;

            row(y(i,j,o,t))=M;

            for l = V_R
              if i!=l
                row(r(i,l,o,t))=M;
              endif
            endfor

            row(s_U(j,i,o,t))=M;
            row(w_U(j,o,t))=-F(o,1)*(t_U(i,j,o)+sigma_u(j,o))/(60*1000);


            constraint_matrix = [constraint_matrix; row];
            RHS = [RHS; 3*M + (((t_U(i,j,o)+sigma_u(j,o))/(60*1000)) * F(o,0)) + ((t_U(i,j,o)+sigma_u(j,o)/(60*1000)) * F(o,1) * w_i(j))];
            Constraint_Types = [Constraint_Types; 'U'];
            U_counter = U_counter + 1;
          endfor
        endfor
        endif
      endfor
    endfor
    disp("41 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%%%%%%%%%%%%%    #constraint 41 B #check
    for i = V_L
      for j = C
        if i!=j
        for o = drone_set
          for t = truck_set
            row = zeros(1,counter);
            row(q(t,j,o))=1;
            row(q(t,i,o))=-1;

            row(y(i,j,o,t))=M;
            for l = V_R
              if i!=l
                row(r(i,l,o,t))=M;
              endif
            endfor

            row(s_U(j,i,o,t))=M;
            row(w_U(j,o,t)) = ((t_U(i,j,o) + sigma_u(j,o)) * F(o,1)) / (60*1000);


            constraint_matrix = [constraint_matrix; row];
            RHS = [RHS; 3*M - ((((t_U(i,j,o)+sigma_u(j,o))/(60*1000)) * F(o,0)) + (((t_U(i,j,o)+sigma_u(j,o))/(60*1000) * F(o,1) * w_i(j))))];
            Constraint_Types = [Constraint_Types; 'U'];
            U_counter = U_counter + 1;
          endfor
        endfor
        endif
      endfor
    endfor
    disp("41 B");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  #constraint 42 A
    for i = C
      for j = C
        for l = V_L
          if i!=j && i!=l && l!=j
          for o = drone_set
             for t = truck_set
                row = zeros(1,counter);

                row(q(t,j,o))=-1;
                row(q(t,i,o))=1;


                row(y(i,j,o,t))=M;
                row(s_U(i,l,o,t))= M;
                row(s_U(j,l,o,t))= M;

                row(w_U(i,o,t))=-(F(o,1)*(t_U(i,j,o)+(sigma_u(j,o))))/(60*1000);

                constraint_matrix = [constraint_matrix; row];
                RHS = [RHS; 3*M +  F(o,0)*(t_U(i,j,o)+(sigma_u(j,o)))/(60*1000)];
                Constraint_Types = [Constraint_Types; 'U'];
                U_counter = U_counter + 1;
             endfor
          endfor
          endif
        endfor
      endfor
    endfor
    disp("42 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%%%%%%%  #constraint 42 B
    for i = C
      for j = C
        for l = V_L
          if i!=j && i!=l && l!=j
          for o = drone_set
             for t = truck_set
                row = zeros(1,counter);

                #q(truck_set,V_L,drone_set)
                row(q(t,j,o))=1;
                row(q(t,i,o))=-1;

                #y(V_L,V_R,drone_set,truck_set)
                row(y(i,j,o,t))=M;

                #s_U(C,V_L,drone_set,truck_set)
                row(s_U(i,l,o,t))= M;
                row(s_U(j,l,o,t))= M;

                #w_U(V,drone_set,truck_set)
                row(w_U(i,o,t))=(F(o,1)*(t_U(i,j,o)+(sigma_u(j,o))))/(60*1000);

                constraint_matrix = [constraint_matrix; row];
                RHS = [RHS; 3*M -  F(o,0)*(t_U(i,j,o)+(sigma_u(j,o)))/(60*1000)];
                Constraint_Types = [Constraint_Types; 'U'];
                U_counter = U_counter + 1;
             endfor
          endfor
          endif
        endfor
      endfor
    endfor
    disp("42 B");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%%%%%
%%   %constraint 43 A
   for i = C
      for j = C
         for o = drone_set
            for t = truck_set
              if i!=j
              row = zeros(1,counter);

              row(y(i,j,o,t))=M;
              row(s_T(j,t))=M;
              row(q(t,i,o)) = -1;

              constraint_matrix = [constraint_matrix; row];
              RHS = [RHS; 2*M - ((F(o,0)*t_U(i,j,o)))/(60*1000)];
              Constraint_Types = [Constraint_Types; 'U'];
              U_counter = U_counter + 1;
              endif
            endfor
         endfor
      endfor
    endfor
    disp("43 A");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);
%%
%%    #constraint 43 b
    for i = C
      for o = drone_set
        for t = truck_set
          row = zeros(1,counter);
          row(q(t,i,o))=-1;
          row(y(i,customers+2,o,t))=M;

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS;M-((F(o,0)*t_U(i,customers+2,o)))/(60*1000)];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
    disp("43 b");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

    #Constraint 46a
    for j = C
       for o = drone_set
          for t = truck_set
             row = zeros(1,counter);
             row(s_T(j,t))=-1;
             for i = V_L
                if i!=j
                   row(y(i,j,o,t))=-1;
                endif
             endfor

             for i = V_L
               if i!=j
                 row(r(i,j,o,t))=M;
               endif
             endfor

             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS;M-2];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
          endfor
       endfor
    endfor
    disp("46a");
    printf("Number of 'S': %d\n", S_counter);
    printf("Number of 'U': %d\n", U_counter);

    #constraint 46b
    for o = drone_set
       for t = truck_set
          row = zeros(1,counter);
          for i = V_L
             row(y(i,customers+2,o,t))=-1;
          endfor

          for i = V_L
             row(r(i,customers+2,o,t))=M;
          endfor

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS;M-1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
       endfor
    endfor

    #constraint 46c
    for i = C
       for o = drone_set
          for t = truck_set
             row = zeros(1,counter);
             row(s_T(i,t))=-1;

             for j = V_R
                if j!=i
                  row(y(i,j,o,t))=-1;
                endif
             endfor

             for j = V_R
                if i!=j
                   row(r(i,j,o,t))=M;
                endif
             endfor

            constraint_matrix = [constraint_matrix; row];
            RHS = [RHS;M-2];
            Constraint_Types = [Constraint_Types; 'U'];
            U_counter = U_counter + 1;
          endfor
       endfor
    endfor

    #constraint 46d
    for o = drone_set
       for t = truck_set
          row = zeros(1,counter);
          for j = V_R
            row(y(1,j,o,t))=-1;
          endfor

          for j = V_R
             row(r(1,j,o,t))=M;
          endfor

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS;M-1];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
       endfor
    endfor

    #fix
    for j = C
      for o = drone_set
         for t = truck_set
            row =  zeros(1,counter);
            for i = V_L
               if i!=j
                  row(r(i,j,o,t))=M;
               endif
            endfor

            for i = V_L
              if i!=j
                row(y(i,j,o,t))=-1;
              endif
            endfor

            row(s_T(j,t))=-1;

           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; M-2];
           Constraint_Types = [Constraint_Types; 'U'];
           U_counter = U_counter + 1;
         endfor
      endfor
    endfor

    #constraint 52
    for i = V_L
      for j = C
         for k = C
            for o = drone_set
               for t = truck_set
                 if i!=j && j!=k && i!=k
                    row = zeros(1,counter);
                    row(p_t(t,j,k))=-1;
                    row(p_t(t,i,k))=M;
                    row(r(i,j,o,t))=M;
                    for l = V_R
                       if l!=k
                          row(r(k,l,o,t))=M;
                       endif
                    endfor
                    constraint_matrix = [constraint_matrix; row];
                    RHS = [RHS; 3*M-1];
                    Constraint_Types = [Constraint_Types; 'U'];
                    U_counter = U_counter + 1;
                 endif
               endfor
            endfor
         endfor
      endfor
    endfor

    #Constraint 53
    for o = drone_set
      for t = truck_set
           row = zeros(1,counter);
           row(y(1,customers+2,o,t))=1;
           constraint_matrix = [constraint_matrix; row];
           RHS = [RHS; 0];
           Constraint_Types = [Constraint_Types; 'S'];
           S_counter = S_counter + 1;
      endfor
    endfor

    #Constraint 54
    for i = C
       for j = C
          if i!=j
            for k = V_L
               if k!=i && k!=j
                  for o = drone_set
                     for t = truck_set
                        row = zeros(1,counter);
                        row(s_U(j,k,o,t))=-1;
                        row(y(i,j,o,t))=M;
                        row(s_U(i,k,o,t))=M;
                        row(s_T(j,t))=-M;

                        constraint_matrix = [constraint_matrix; row];
                        RHS = [RHS; 2*M-1];
                        Constraint_Types = [Constraint_Types; 'U'];
                        U_counter = U_counter + 1;
                     endfor
                  endfor
               endif
            endfor
          endif
       endfor
    endfor

    #constraint 55
    for i = C
       for o = drone_set
          for t = truck_set
             row = zeros(1,counter);
             row(t_UA(i,o,t))=1;
             row(t_UA(customers+2,o,t))=-1;
             for j = V_L
               row(r(j,customers+2,o,t))=M;
             endfor

             constraint_matrix = [constraint_matrix; row];
             RHS = [RHS; M];
             Constraint_Types = [Constraint_Types; 'U'];
             U_counter = U_counter + 1;
          endfor
       endfor
    endfor




    #constraint that prevents all energy capacities from being greater than the max batter capacity
    for o = drone_set
      for t = truck_set
        for i = V_L
          row = zeros(1,counter);
          row(q(t,i,o))=1;

          constraint_matrix = [constraint_matrix; row];
          RHS = [RHS; e(o)];
          Constraint_Types = [Constraint_Types; 'U'];
          U_counter = U_counter + 1;
        endfor
      endfor
    endfor
  disp("THE TOTAL NUMBER OF CONSTRAINTS");
  disp(U_counter+S_counter);

  lb = zeros(1, big_T);
  ub = Inf * ones(1, big_T);

  for i = 1:big_T
    if variable_types(i) == 'B'
        variable_types(i) = 'I';
        ub(i) = 1;
    elseif variable_types(i) == 'I'
        ub(i) = customers+2;
    elseif variable_types(i) == 'C'
        ub(i) = 10000000;
    endif
  endfor


  sense = 1;

  ##preallocate of size constraint matrix
  num_U = sum(Constraint_Types == 'U');
  num_S = sum(Constraint_Types == 'S');

  ineq_matrix = zeros(num_U, size(constraint_matrix, 2));
  eq_matrix   = zeros(num_S, size(constraint_matrix, 2));

  ineq_count = 1;
  ineq_RHS = [];

  eq_count = 1;
  eq_RHS = [];

  for i = 1:size(constraint_matrix, 1)
     if Constraint_Types(i) == 'U'
        ineq_matrix(ineq_count, :) = constraint_matrix(i, :);
        ineq_RHS = [ineq_RHS; RHS(i)];
        ineq_count = ineq_count + 1;
     elseif Constraint_Types(i) == 'S'
        eq_matrix(eq_count, :) = constraint_matrix(i, :);
        eq_RHS = [eq_RHS; RHS(i)];
        eq_count = eq_count + 1;
    end
  end
  disp(size(ineq_matrix));
  disp(size(eq_matrix));
  [solution_vector, fval]=cplex_solve(1,c,ineq_matrix,ineq_RHS,eq_matrix, eq_RHS,lb,ub,variable_types,'LP_truck_drone',cd(),0,1);


  fprintf("Objective Value %f\n",fval);
  interpret_solution(solution_vector,customers, truck_set,drone_set, V_L, V_R, C, x, s_U, r, y, t_TA,t_TL,t_UA,t_UL,w_i,Q,q,e,s_T,w_U);

  figure;
  plot_customer_points(customer_points, solution_vector, ...
                              x, y, V_L, V_R, truck_set, drone_set);

endfunction


function [C, C_T, C_U, C_F] = customer_sets(customers)
  C = 2:customers+1;
  proportions = [0.33, 0.33, 0.34];

  % Step 1: Randomly permute the customer vector
  C_perm = C(randperm(length(C)));

  % Step 2: Determine how many customers go in each group
  n = length(C);
  counts = round(proportions * n);

  % Ensure total count equals n by correcting rounding errors
  diff = n - sum(counts);
  [~, idx] = max(proportions);  % assign difference to the largest group
  counts(idx) = counts(idx) + diff;

  % Step 3: Assign customers
  C_T = C_perm(1:counts(1));
  C_U = C_perm(counts(1)+1:counts(1)+counts(2));
  C_F = C_perm(counts(1)+counts(2)+1:end);


endfunction



function [customer_points, distance_matrix] = make_distance_matrix(num_customers, width, height)

    total_nodes   = num_customers + 2;
    node_coords   = zeros(total_nodes, 2);

    node_coords(1, :)   = [0, 0];
    node_coords(end, :) = [0, 0];



    x = (2*rand(num_customers,1) - 1) * width;
    node_coords(2 : total_nodes-1, 1) = x;

    y = (2*rand(num_customers,1) - 1) * height;
    node_coords(2 : total_nodes-1, 2) = y;


    customer_points = node_coords;

    distance_matrix = zeros(total_nodes);
    for i = 1 : total_nodes
        for j = 1 : total_nodes
            dx = node_coords(i,1) - node_coords(j,1);
            dy = node_coords(i,2) - node_coords(j,2);
            distance_matrix(i,j) = hypot(dx, dy);
        end
    end
end



function plot_customer_points(customer_points, solution_vector, ...
                              x, y_idx, V_L, V_R, truck_set, drone_set)

    if size(customer_points,2) ~= 2
        error('customer_points must be an n‑by‑2 matrix of [x y] coordinates.');
    end

    clf; hold on;

    % Plot customers (blue) and depots (red)
    plot(customer_points(2:end-1,1), customer_points(2:end-1,2), ...
         'bo', 'MarkerFaceColor','b', 'MarkerSize',6);
    plot(customer_points([1,end],1), customer_points([1,end],2), ...
         'rs', 'MarkerFaceColor','r', 'MarkerSize',8);

    % Label all nodes with their index
    for k = 1:size(customer_points,1)
        text(customer_points(k,1)+0.5, customer_points(k,2)+0.5, ...
             num2str(k), 'FontSize',10, 'Color','k');
    end

    % === Truck arcs (solid black arrows) =======================================
    for t = truck_set
        for i = V_L
            for j = V_R
                if i ~= j && solution_vector(x(i,j,t)) > 0.8
                    p1 = customer_points(i,:);
                    p2 = customer_points(j,:);
                    dx = p2(1) - p1(1);
                    dy = p2(2) - p1(2);
                    quiver(p1(1), p1(2), dx, dy, 0, ...
                           'Color', 'k', 'LineWidth', 2, ...
                           'MaxHeadSize', 0.05, 'AutoScale', 'off');
                end
            end
        end
    end

    % === Drone arcs (solid green arrows) =====================================
    for o = drone_set
        for t = truck_set
            for i = V_L
                for j = V_R
                    if i ~= j && solution_vector(y_idx(i,j,o,t)) > 0.8
                        p1 = customer_points(i,:);
                        p2 = customer_points(j,:);
                        dx = p2(1) - p1(1);
                        dy = p2(2) - p1(2);
                        quiver(p1(1), p1(2), dx, dy, 0, ...
                               'Color', 'g', 'LineWidth', 1.5, ...
                               'MaxHeadSize', 0.03, 'AutoScale', 'off');
                    end
                end
            end
        end
    end

    title('Truck (black) and Drone (green) Paths');
    xlabel('X coordinate');  ylabel('Y coordinate');
    axis equal; grid on; hold off;
end




function [w_i] = Generate_Weight(C,C_T,C_U,C_F,Q,Z,num_truck,num_drone)
  weight_ratio = 0.15;
  max_id = max(C);
  w_i = zeros(1, max_id);

  truck_max = weight_ratio * min(Z(1:num_truck));
  drone_max =  weight_ratio * max(Q(1:num_drone));
  flex_max = weight_ratio * min(min(Q),min(Z));

  disp(C);
  for i = C
    if ismember(i, C_T)
      w_i(i)=rand()*truck_max;
    endif

    if ismember(i, C_U)
      disp(i);
      t=rand()*drone_max;
      w_i(i)=t;
      disp(t);
    endif

    if ismember(i,C_F)
      w_i(i)=rand()*flex_max;
    endif
  endfor

endfunction





function [a] = F(drone_type,coeff)
   if drone_type == 1
      if coeff==0
        a = 783.3; #a_0
      else
        a = 354.7849; #a_1
      endif
   endif
   if drone_type == 2
     if coeff==0
       a = 1093.6; #a_0
     endif
     if coeff==1
       a = 232.3; #a_1
     endif
   endif
endfunction




function [] = interpret_solution(solution_vector,customers, truck_set, drone_set, V_L, V_R, C, x, s_U, r, y, t_TA, t_TL, t_UA, t_UL,w_i,Q,q,e,s_T,w_U)
  for t = truck_set
    fprintf("Truck %d route :\n",t)
    for i = V_L
      for j = V_R
        if i != j && solution_vector(x(i, j, t)) > 0.8
           arrival = solution_vector(t_TA(j,t));
           depature = solution_vector(t_TL(i,t));
           fprintf("[%d(%f) %d(%f]\n",i,depature,j,arrival)
        endif
      endfor
    endfor
  endfor


  % Print drone routes
  for o = drone_set
    for t = truck_set
      fprintf('Drone %d  |  Truck %d:\n', o, t);
      fprintf('Capacity(%f kg) \n', Q(o));
      disp("services customers:");
      for i = C
        for j = V_L
          if i!=j
            if solution_vector(s_U(i,j,o,t))>0.8
              fprintf("Launches from node %d to service customer %d ", j, i);
              fprintf("parcel weight (%f) ", w_i(i));
              for temp = V_R
                if j!=temp && solution_vector(r(j,temp,o,t))>0.8
                  fprintf("retreived at node %d\n",temp);
                endif
              endfor
            endif
          endif
        endfor
      endfor


      % find and print every edge for this (o,t) pair
      for i = V_L
        for j = V_R
          if i ~= j && solution_vector(y(i, j, o, t)) > 0.8
            #time
            arrival_time = solution_vector(t_UA(j,o,t));
            departure_time = solution_vector(t_UL(i,o,t));

            #weight
            departure_weight = solution_vector(w_U(i,o,t));
            arrival_weight = solution_vector(w_U(j,o,t));

            #energy
            depature_energy = solution_vector(q(t,i,o));
            if j!=customers + 2
               arrival_energy = solution_vector(q(t,j,o));
            else
               arrival_energy = 0;
            endif


            fprintf('  [%d(depature time: %f weight: %f energy: %f), %d(arrival time: %f weight: %f energy: %f)]\n', i,departure_time,departure_weight,depature_energy, j,arrival_time,arrival_weight,arrival_energy);
         end
        end
      end
      fprintf('\n');
    end
  end


endfunction



