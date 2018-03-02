function y = mismatch_locations(N_genome_i)
%make sure the data is formatted correctly
flip = size(N_genome_i);
if flip(2)>flip(1)
    N_genome_i = N_genome_i';
end
%convert to 0 1
% Define 0=G,1=C,2=A,3=T
t = size(N_genome_i);
for i = 1:max(t)
   %match
   if (N_genome_i(i,1) == 0 && N_genome_i(i,2) == 1) || (N_genome_i(i,2) == 0 && N_genome_i(i,1) == 1)
       newG(i,1) = 0;
       newG(i,2) = 0;
   elseif (N_genome_i(i,1) == 2 && N_genome_i(i,2) == 3) || (N_genome_i(i,2) == 2 && N_genome_i(i,1) == 3)
       newG(i,1) = 1;
       newG(i,2) = 1;
   %mismatch
   elseif (N_genome_i(i,1) == 1 && N_genome_i(i,2) == 3) || (N_genome_i(i,2) == 1 && N_genome_i(i,1) == 3)
       newG(i,1) = 1;
       newG(i,2) = 0;
   elseif (N_genome_i(i,1) == 0 && N_genome_i(i,2) == 3) || (N_genome_i(i,2) == 0 && N_genome_i(i,1) == 3)
       newG(i,1) = 1;
       newG(i,2) = 0;
   elseif (N_genome_i(i,1) == 1 && N_genome_i(i,2) == 2) || (N_genome_i(i,2) == 1 && N_genome_i(i,1) == 2)
       newG(i,1) = 0;
       newG(i,2) = 1;
   else
       newG(i,1) = 0;
       newG(i,2) = 1;
   end
end
N_genome_i = newG;
complimentary_strand = N_genome_i(:,1);
non_complimentary_strand = N_genome_i(:,2);
mis_locations = zeros(size(N_genome_i(:,1)));
ii = 1;
while ii <= length(complimentary_strand)
    if complimentary_strand(ii) == non_complimentary_strand(ii)
        mis_locations(ii) = 1; %match
    else
        mis_locations(ii) = 0; %mismatch
    end
ii = ii+1;    
end
y = mis_locations;
end