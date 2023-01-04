class reaclib_entry:
     """ class for individual reaclib ENTRIES, containing reaction
         information and the 7 fit coefficients. 
         The reactands are stored as nucleus objects.
         Logic Variables  resonant, reverse and weak corresponding 
         to the respective REACLIB markers.
         Reaclib_entries for the same reactions are to be summarized to 
         reaclib_reactions.
     """
     a=[0.0]*7     
     reactands=5
     reactands=["====="]*6
     label="===="
     weak=False
     resonant=False
     reverse=False
     qval=0.0
     chapter=0
     def __init__(self,chapter,a,nucs,label,weak,resonant,reverse,qval):
         """ -chapher is the reaclib chapter (integer), e.g. 4 for (n,g), 2 for (g,n)
         -a is list/array with 7 entries
         -nucs is a list of involved nuclei with 6 entries, p/n/he4 always first
         e.g ga84(n,g)ga85 : ["n","ga84","ga85","=====","=====","====="
         -label is a four character indicator for the source of the reaction
         -weak : True or False for weak reaction
         -resonant : True of False for resonant reaction
         -reverse : is it a derived reverse reaction (do we need to include partition functions?) True/False
         -qval : Q-value (MeV)  
          ]
         """      
         for j in range(0,6):
           self.reactands[j]="====="
         self.a=a
         self.chapter=chapter
         for j in range(0,len(nucs)):
           self.reactands[j]=nucs[j]
         
# guarantees that there are always 6 entries
         self.label=label
         self.weak=weak
         self.resonant=resonant
         self.reverse=reverse  
         self.qval=qval
     def get_header(self):
         line="     "
         for j in range(0,len(self.reactands)):
            reactand=self.reactands[j]
            line=line+reactand.replace("="," ").rjust(5)
         line=line+"        "
         line=line+self.label
         if self.resonant:
             line=line+"r"
         else:
             line=line+" "
         if self.reverse:
             line=line+"v"
         else:
              line=line+" "
         line=line+"   "
         line=line+str("%12.5e" % self.qval)        
         line=line+"          \n"
         return(line)
     def print_coeffs(self):
         line1=""
         line2=""
         for j in range(0,4):
             line1=line1+str("%13.6e" % self.a[j])
         line1=line1+"                      \n"
         for j in range(4,7):
             line2=line2+str("%13.6e" % self.a[j])
         line2=line2+"                                   \n"                 
         return(line1,line2)
     def print_to_file(self,filename):
         header=self.get_header()
         body=self.print_coeffs()
         
         testout=open(filename,"a")
         testout.write(header)
         testout.write(body[0])
         testout.write(body[1])
         
