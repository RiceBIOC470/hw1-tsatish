% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
x = '3'; y= '5'; %strings
x = 3; y = '5'; %mixed

%your code goes here

if ~isnumeric(x)
    x = str2num(x);
end
if ~isnumeric(y)
    y = str2num(y);
end
result = x + y;

%output your answer
disp(result)

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 

N = 500; % define sequence length
random_data = randi(4,1,N);
rand_seq = char(zeros(1,N));
for i = 1:N
    switch random_data(1,i)
        case 1
            rand_seq(1,i) = 'A';
        case 2
            rand_seq(1,i) = 'T';
        case 3
            rand_seq(1,i) = 'G';
        case 4
            rand_seq(1,i) = 'C';
    end
end

%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

%Find first index of all possible start codons
start_locations = strfind(rand_seq, 'ATG');
%Find first index of all possible stop codons
end_locations = cat(2, strfind(rand_seq, 'TAA'), strfind(rand_seq, 'TGA'), strfind(rand_seq, 'TAG'));
ORF_lengths = [];
%For each start codon find the next stop codon and use it to calculate ORF
%length
for x = 1:length(start_locations)
    potential_stop_points = [];
    for y = 1:length(end_locations)
        if end_locations(y) > start_locations(x)
            potential_stop_points(end+1) = end_locations(y);
        end
    end
    %Select stop codon immediately after selected start codon
    actual_stop_point = min(potential_stop_points);
    %Getting the ORF length:
    %Add 2 because the strfind gets the location of the first letter, so we
    %would miss the last two bases of the stop codon
    this_ORF_length = actual_stop_point - start_locations(x) + 2;
    ORF_lengths(end+1) = this_ORF_length;
end

disp(max(ORF_lengths));
    
%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
count = 0;
for temp = 1:1000
    N = 500; % define sequence length
    random_data = randi(4,1,N);
    rand_seq = char(zeros(1,N));
    for i = 1:N
        switch random_data(1,i)
            case 1
                rand_seq(1,i) = 'A';
            case 2
                rand_seq(1,i) = 'T';
            case 3
                rand_seq(1,i) = 'G';
            case 4
                rand_seq(1,i) = 'C';
        end
    end
   
    %Find first index of all possible start codons
    start_locations = strfind(rand_seq, 'ATG');
    %Find first index of all possible stop codons
    end_locations = cat(2, strfind(rand_seq, 'TAA'), strfind(rand_seq, 'TGA'), strfind(rand_seq, 'TAG'));
    ORF_lengths = [];
    %For each start codon find the next stop codon and use it to calculate ORF
    %length
    for x = 1:length(start_locations)
        potential_stop_points = [];
        for y = 1:length(end_locations)
            if end_locations(y) > start_locations(x)
                potential_stop_points(end+1) = end_locations(y);
            end
        end
        %Select stop codon immediately after selected start codon
        if ~isempty(potential_stop_points)
            actual_stop_point = min(potential_stop_points);
            %Getting the ORF length:
            %Add 2 because the strfind gets the location of the first letter, so we
            %would miss the last two bases of the stop codon
            this_ORF_length = actual_stop_point - start_locations(x) + 2;
            ORF_lengths(end+1) = this_ORF_length;
        end        
    end
    
    if max(ORF_lengths) > 50
        count = count+1;
    end
end

disp(count/1000);


%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
prob_of_50 = [];
max_length = 1000;
for current_N = 1:max_length
    count = 0;
    sim_length=1000;
    for temp = 1:sim_length
        N = current_N; % define sequence length
        random_data = randi(4,1,N);
        rand_seq = char(zeros(1,N));
        for i = 1:N
            switch random_data(1,i)
                case 1
                    rand_seq(1,i) = 'A';
                case 2
                    rand_seq(1,i) = 'T';
                case 3
                    rand_seq(1,i) = 'G';
                case 4
                    rand_seq(1,i) = 'C';
            end
        end

        %Find first index of all possible start codons
        start_locations = strfind(rand_seq, 'ATG');
        %Find first index of all possible stop codons
        end_locations = cat(2, strfind(rand_seq, 'TAA'), strfind(rand_seq, 'TGA'), strfind(rand_seq, 'TAG'));
        ORF_lengths = [];
        %For each start codon find the next stop codon and use it to calculate ORF
        %length
        for x = 1:length(start_locations)
            potential_stop_points = [];
            for y = 1:length(end_locations)
                if end_locations(y) > start_locations(x)
                    potential_stop_points(end+1) = end_locations(y);
                end
            end
            %Select stop codon immediately after selected start codon
            if ~isempty(potential_stop_points)
                actual_stop_point = min(potential_stop_points);
                %Getting the ORF length:
                %Add 2 because the strfind gets the location of the first letter, so we
                %would miss the last two bases of the stop codon
                this_ORF_length = actual_stop_point - start_locations(x) + 2;
                ORF_lengths(end+1) = this_ORF_length;
            end        
        end

        if max(ORF_lengths) > 50
            count = count+1;
        end
    end
    prob_of_50(current_N) = count/sim_length;
end

plot(1:max_length, prob_of_50)

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

% When N is below 50 the probability should be 0, and when N is just above
% 50 the probability should be vanishingly small. Only when N becomes
% extremely large should the probability start becoming significant.
% Intuitively it seems the increase in the probability should be gradual
% and not really linear. The graph seems to match all these
% characteristics.

%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

filename = 'qPCRdata.txt';
fileID = fopen(filename,'r');
formatSpec = '%*s%*s%s%*s%f%*s%*s%*s%[^\n]';
% Need rows 3-74, 72 lines to read, 2 header lines, tab delimited
data = textscan(fileID, formatSpec, 72, "Delimiter", "\t","HeaderLines", 2);
Cp_vector = cell2mat(data(1,2));

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 

plate = zeros(6,12);
count = 0;
for y = 1:6
    for x= 1:12
        count = count+1;
        plate(y,x) = Cp_vector(count);
    end
end

% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 

%Questions: are you supposed to be averaging the Cp values as I am? Are you
%supposed to display the normalization gene and condition 1 on the map?

heat_map = zeros(6,4);
CpN0 = (plate(1,10) + plate(1,11) + plate(1,12))/3;
y_count = 1;
for y = 1:3:11
    Cp0 = (plate(1,y) + plate(1,y+1) + plate(1,y+2))/3;
    for x = 1:6
        CpX = (plate(x,y) + plate(x,y+1) + plate(x,y+2))/3;
        CpNX = (plate(x,10) + plate(x,11) + plate(x,12))/3;
        fold_change = 2^(Cp0 - CpX - (CpN0 - CpNX));
        heat_map(x,y_count) = fold_change;
    end
    y_count = y_count+1;
end
        
heatmap(heat_map);


%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


