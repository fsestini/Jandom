//! @tags complexity
//! @citations AliasDFG10

//Generated by c2fsm -cut -nodiv -int 
model ax {
var i,j,n;
//parameters n;
states stop,start,lbl_7_1,cut;
transition t_24 :={
  from  := start;
  to    := lbl_7_1;
  guard := (2 <= n);
  action:= i' = 0, j' = 1;
};
transition t_26 :={
  from  := start;
  to    := stop;
  guard := (n <= 1);
  action:= i' = 1, j' = 0;
};
transition t_18 :={
  from  := lbl_7_1;
  to    := lbl_7_1;
  guard := (j+2 <= n);
  action:= j' = j+1;
};
transition t_19 :={
  from  := lbl_7_1;
  to    := cut;
  guard := ( (i+3 <= n) && (n <= j+1) );
  action:= i' = i+1;
};
transition t_20 :={
  from  := lbl_7_1;
  to    := stop;
  guard := ( (n <= i+2) && (n <= j+1) );
  action:= i' = i+1;
};
transition t_21 :={
  from  := cut;
  to    := lbl_7_1;
  guard := (2 <= n);
  action:= j' = 1;
};
transition t_22 :={
  from  := cut;
  to    := cut;
  guard := ( (i+3 <= n) && (n <= 1) );
  action:= i' = i+1, j' = 0;
};
transition t_23 :={
  from  := cut;
  to    := stop;
  guard := ( (n <= i+2) && (n <= 1) );
  action:= i' = i+1, j' = 0;
};
}
strategy dumb {
    Region init := { state = start };
}

