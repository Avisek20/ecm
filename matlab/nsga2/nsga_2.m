function nsga_2(pop,gen)
%% function nsga_2(pop,gen)
% is a multi-objective optimization function where the input arguments are 
% pop - Population size
% gen - Total number of generations
% 
% This functions is based on evolutionary algorithm for finding the optimal
% solution for multiple objective i.e. pareto front for the objectives. 
% Initially enter only the population size and the stoping criteria or
% the total number of generations after which the algorithm will
% automatically stopped. 
%
% You will be asked to enter the number of objective functions, the number
% of decision variables and the range space for the decision variables.
% Also you will have to define your own objective funciton by editing the
% evaluate_objective() function. A sample objective function is described
% in evaluate_objective.m. Kindly make sure that the objective function
% which you define match the number of objectives that you have entered as
% well as the number of decision variables that you have entered. The
% decision variable space is continuous for this function, but the
% objective space may or may not be continuous.
%
% Original algorithm NSGA-II was developed by researchers in Kanpur Genetic
% Algorithm Labarotary and kindly visit their website for more information
% http://www.iitk.ac.in/kangal/


%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

%% Simple error checking
% Number of Arguments
% Check for the number of arguments. The two input arguments are necessary
% to run this function.
if nargin < 2
    error('NSGA-II: Please enter the population size and number of generations as input arguments.');
end
% Both the input arguments need to of integer data type
if isnumeric(pop) == 0 || isnumeric(gen) == 0
    error('Both input arguments pop and gen should be integer datatype');
end
% Minimum population size has to be 20 individuals
if pop < 20
    error('Minimum population for running this function is 20');
end
if gen < 5
    error('Minimum number of generations is 5');
end
% Make sure pop and gen are integers
pop = round(pop);
gen = round(gen);
%% Objective Function
% The objective function description contains information about the
% objective function. M is the dimension of the objective space, V is the
% dimension of decision variable space, min_range and max_range are the
% range for the variables in the decision variable space. User has to
% define the objective functions using the decision variables. Make sure to
% edit the function 'evaluate_objective' to suit your needs.
[M, V, min_range, max_range] = objective_description_function();

%% Initialize the population
% Population is initialized with random values which are within the
% specified range. Each chromosome consists of the decision variables. Also
% the value of the objective functions, rank and crowding distance
% information is also added to the chromosome vector but only the elements
% of the vector which has the decision variables are operated upon to
% perform the genetic operations like corssover and mutation.
global chromosome;
chromosome = initialize_variables(pop, M, V, min_range, max_range);

% cen1= [0.347741462581702,0.283567427832546,0.563040131460230,0.446408994344749,-0.961960710258747,-1,0.994487586200187,-0.945717989669797;-0.0168622307170359,-0.0868140170984634,0.111968110640444,0.0846241919914952,-0.161359594808892,-0.199505491299369,0.238884703933920,0.300331267074125;-0.0725806121660444,-0.0303719159695115,0.0393082889816514,0.0580938786335623,-0.0452308977359758,-0.0398367869236616,0.261781569536408,0.295117636468225;-0.0514534650907537,-0.104199834848636,0.115797982901305,0.0942524884239095,-0.0227128645842385,-0.118298614577129,0.353929890863459,0.375704205023120;0.279749786514881,0.146540391455874,0.449479344332190,0.318991960825518,-0.906605211602913,-0.974210305418633,1,-0.942305419003487;0.345466100758285,0.0780811461977184,0.458157175764478,0.400260265689594,-0.942117212806172,-0.990012307995070,1,-0.884282517679079;0.271259937633361,0.103353774949189,0.159405358175828,0.184768117605406,-0.823050756041641,-0.904154294140568,1,-0.892250065701488;0.201170059704918,0.0900059895236135,0.360749326142351,0.270513432890786,-0.873556940333273,-0.923927991963644,1.00000000000000,-0.883367829631880;-0.0564040204947056,-0.0203036749061046,0.0570501931101689,0.0612424948161064,-0.135063503566461,-0.0242118150349526,0.207026843149427,0.416447344882664;0.338879658738897,0.289210926253233,0.537768168746553,0.400583259068421,-1,-0.983690331184886,1,-0.953278020490789;0.219382444217144,0.221548909319488,-0.102933001604380,-0.000513386374182842,-0.982218328794181,1,-1,-1;0.0795266764424534,-0.0150813541798465,0.0155691074168549,0.106236533207931,0.140004546910722,0.150537424304542,0.391132941462564,0.267846501167297;-0.0634034329587275,-0.0254186539666016,0.0451020491157070,0.0600284584307461,-0.0103111527882292,-0.0391291152109890,0.261729710534134,0.261223434177713;0.292456659684256,0.117188095493739,0.172206387629743,0.181807111297950,-0.863355929043478,-0.892376379851157,1,-0.891042680605868;-0.0271450487590152,-0.0907193945060832,0.108257523313161,0.0992161232995318,-0.0321608797257923,-0.116789989996397,0.235621382474396,0.319592945344050;-0.0668389260206328,-0.0249368308613040,0.0391279574918094,0.0569036073024049,-0.0657507838755830,-0.0295644606945320,0.269988412547101,0.345657285663382;-0.0539523514653436,-0.0214504334660086,0.0603378130084834,0.0648485275062664,-0.0927287945912506,0.0252980719469418,0.274371692363932,0.374900542921417;-0.0514109343574368,-0.0969845710212997,0.121708749932982,0.0968577497247347,-0.0501318952225772,-0.117974418991746,0.239282898372002,0.409646920182875;0.332343381782613,0.272190741649253,0.436698951640745,0.398326218330094,-1,-0.996429150431087,1,-0.925699470788662;0.199273078759075,0.120668687400767,0.369787043765832,0.260209311649464,-0.905546128205828,-0.953317264050660,1,-0.879915156205609;0.329822997994394,0.274762334378010,0.467700116719999,0.397481115531735,-0.962062701039705,-0.995656444149519,1,-0.961773040041615;-0.102334336144757,0.00500835999743859,0.0523639884589982,0.00745736080452205,0.0947713501480307,0.0644800941500912,0.364147414975531,0.202421110572315;-0.0793499713534637,-0.0151419834877965,0.0244031585492770,0.0555592684395981,-0.115863542117295,-0.0541651225150282,0.268341714709220,0.291608666392696;-0.0300618967666558,-0.0829062771566220,0.117704717608076,0.0913231104844552,-0.0709095646981957,-0.110978526437562,0.240103051042311,0.334780722023526;0.306105471270879,0.251371545288651,0.444127726897502,0.401300152525146,-0.977185344577861,-0.980577079259779,1,-0.878886646523766;0.291466098706427,0.138319575708879,0.346706563426243,0.329948762428455,-0.924088615150134,-1,1,-0.903872658062760;0.235325517013562,0.167929705786681,0.347861812429428,0.321806062769807,-0.979042824021038,-0.999872703885157,0.997146270975129,-0.915402997206483;-0.0403764023912594,0.223790986907836,-0.0117236115375027,-0.0823636652988713,-0.915832821108414,1,-0.850297097328357,-0.841528123865006;-0.0269977798105707,-0.0849927296634400,0.115015905803329,0.0875768789251554,-0.0324821470668795,-0.151959589324726,0.243418626494895,0.313104287811519;-0.0289160564575581,-0.0873555456366672,0.116115414921658,0.0863175267072469,-0.0564309553524050,-0.135388123005686,0.237249322017440,0.323663993797944];
% 
% chromosome= cen1;

%% Sort the initialized population
% Sort the population using non-domination-sort. This returns two columns
% for each individual which are the rank and the crowding distance
% corresponding to their position in the front they belong. At this stage
% the rank and the crowding distance for each chromosome is added to the
% chromosome vector for easy of computation.
chromosome = non_domination_sort_mod(chromosome, M, V);

%% Start the evolution process
% The following are performed in each generation
% * Select the parents which are fit for reproduction
% * Perfrom crossover and Mutation operator on the selected parents
% * Perform Selection from the parents and the offsprings
% * Replace the unfit individuals with the fit individuals to maintain a
%   constant population size.

for i = 1 : gen
    
    % Select the parents
    % Parents are selected for reproduction to generate offspring. The
    % original NSGA-II uses a binary tournament selection based on the
    % crowded-comparision operator. The arguments are 
    % pool - size of the mating pool. It is common to have this to be half the
    %        population size.
    % tour - Tournament size. Original NSGA-II uses a binary tournament
    %        selection, but to see the effect of tournament size this is kept
    %        arbitary, to be choosen by the user.
    pool = round(pop/2);
    tour = 2;
    % Selection process
    % A binary tournament selection is employed in NSGA-II. In a binary
    % tournament selection process two individuals are selected at random
    % and their fitness is compared. The individual with better fitness is
    % selcted as a parent. Tournament selection is carried out until the
    % pool size is filled. Basically a pool size is the number of parents
    % to be selected. The input arguments to the function
    % tournament_selection are chromosome, pool, tour. The function uses
    % only the information from last two elements in the chromosome vector.
    % The last element has the crowding distance information while the
    % penultimate element has the rank information. Selection is based on
    % rank and if individuals with same rank are encountered, crowding
    % distance is compared. A lower rank and higher crowding distance is
    % the selection criteria.
    parent_chromosome = tournament_selection(chromosome, pool, tour);

    % Perfrom crossover and Mutation operator
    % The original NSGA-II algorithm uses Simulated Binary Crossover (SBX) and
    % Polynomial  mutation. Crossover probability pc = 0.9 and mutation
    % probability is pm = 1/n, where n is the number of decision variables.
    % Both real-coded GA and binary-coded GA are implemented in the original
    % algorithm, while in this program only the real-coded GA is considered.
    % The distribution indeices for crossover and mutation operators as mu = 20
    % and mum = 20 respectively.
    mu = 20;
    mum = 20;
    offspring_chromosome = ...
        genetic_operator(parent_chromosome, ...
        M, V, mu, mum, min_range, max_range);

    % Intermediate population
    % Intermediate population is the combined population of parents and
    % offsprings of the current generation. The population size is two
    % times the initial population.
    
    [main_pop,temp] = size(chromosome);
    [offspring_pop,temp] = size(offspring_chromosome);
    % temp is a dummy variable.
    clear temp
    % intermediate_chromosome is a concatenation of current population and
    % the offspring population.
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = ...
        offspring_chromosome;

    % Non-domination-sort of intermediate population
    % The intermediate population is sorted again based on non-domination sort
    % before the replacement operator is performed on the intermediate
    % population.
    intermediate_chromosome = ...
        non_domination_sort_mod(intermediate_chromosome, M, V);
    % Perform Selection
    % Once the intermediate population is sorted only the best solution is
    % selected based on it rank and crowding distance. Each front is filled in
    % ascending order until the addition of population size is reached. The
    % last front is included in the population based on the individuals with
    % least crowding distance
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    %if ~mod(i,100)
        %clc
        %fprintf('%d generations completed\n',i);
    %end
end

%% Result
% Save the result in ASCII text format.
%save solution.txt chromosome -ASCII
assignin('base','chromosome',chromosome);
%% Visualize
% The following is used to visualize the result if objective space
% dimension is visualizable.
%if M == 2
%    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
%elseif M ==3
%    plot3(chromosome(:,V + 1),chromosome(:,V + 2),chromosome(:,V + 3),'*');
%end
end   
