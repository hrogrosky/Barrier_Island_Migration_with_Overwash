function [Hstar]=AeolianTransport03312021(Hstar,delta,L,W,Windspeed,Windmin,StormThreshold,WindDir,PC)

%Aeolian Transport is triggered by sufficient wind speeds chosen in the MainCode.
%AT determines the number of cells that the sediment will move based on the
%current windspeed:
%                           windspeed < 6 move 1 cell
%                           6<= windspeed < 11 move 2 cells
%                           11<= windspeed < 16 move 3 cells

%When each cell is polled the plant density on that cell is accounted for,
%but it is assumed that after initial saltation the sediment is moving
%freely enough so as not to require polling subsequent cells' plant cover.
%The angle between the current cell and the next potential cell is checked
%for every step 1, 2, and 3.

%each cell has a probability of moving either in the same direction as the
%WindDirection input or one of the two off-directions 
% (i.e. if WindDir is North: 50% move North, 25% move NorthEast, 25% move NorthWest


MvF=.8; %movement probability factor (turn up for more movement, down for less) 

Hp=padarray(Hstar(:,:,1),[3 3]);    %padding Hstar with three estra rows/cols in each direction - avoids index errors
                                        %cells moved into padded region are
                                        %deleted when Hstar is redeclared
                                        %after each loop
[n1 n2]=size(Hstar(:,:,1));         %array indexes for loop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('I am running AT! ')


%can move up to three cells based on windspeed (distance not amount)
MvCnt=1*(Windspeed>=6)+1*(Windspeed>=9)+1*(Windspeed>=12);   

%***about MarshFlag*** - we are currently allowing the wind to affect cells in
%the swamp, but previously we were flagging them so that they would not be
%moved. If you do not want the wind to move swamp cells then change the
%rows labeled below

for i=4:n1+3
    Rnow=Hp(i,:);
    for j=n2+3:-1:4
        MarshFlag=0;        %flagging cells in the marsh
        if Rnow(j)>=0
            if (Rnow(j)> 0)&& Rnow(j)<=1
                if j==1
                    if Rnow(j)<Rnow(j+1) && Rnow(j+1)<Rnow(j+2)
                        MarshFlag=0;                                %<---change to MarshFlag=1; if wish to ignore marsh cells
                    end
                elseif j==n2
                    if Rnow(j)>Rnow(j-1) && Rnow(j-1)>Rnow(j-2)
                        MarshFlag=0;                                %<---change to MarshFlag=1; if wish to ignore marsh cells
                    end
                else
                    if Rnow(j)<Rnow(j+1) && Rnow(j-1)<Rnow(j)
                        MarshFlag=0;                                %<---change to MarshFlag=1; if wish to ignore marsh cells
                    end    
                end
            end
            
            if Rnow(j)>0 && MarshFlag==0        
                    %probability for RNG use
                    P=0.5;      %P splits the likelihood of moving in primary wind direction or the two alternative directions
                    Rand_dir=rand;  %for probability check
                    MoveChance=(1-PC(i,j))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;   %cells chance of moving, factors in plants and windpseed ranges
                    dist=round(MvCnt*MoveChance,0); %number of cells a slab can potentially move

                    %%%%%%%%%%%%%%%%%%%%%%%
                    %% DEPOSITION RULES: %%
                    %%%%%%%%%%%%%%%%%%%%%%%

                    %Wind Direction
                %1=N  2=NE  3=E  4=SE  5=S  6=SW  7=W  8=NW 

                    Nhd=Hp(i-3:i+3,j-3:j+3); %7x7 moving window for movement up to three cells
                    PCNhd=PC(i-3:i+3,j-3:j+3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %NORTH BLOCK
                    
                    %NORTH - NW
                    
                    if WindDir==1 && rand<MoveChance       
                        if Rand_dir<((1-P)/2) %&& (i>3) && (j>3) %move to the port slab I don't think I need to worry about these if I use the padded array
                            DirV=[Nhd(4,4), Nhd(3,3), Nhd(2,2), Nhd(1,1)];  %vector of neighborhood Hstar values in the direction of chosen movement
                            
            %next line for re-calculating chance of moving from new cell - if we wish to poll for plant cover at each step at some later time
                            PCNE=[PCNhd(4,4), PCNhd(3,3), PCNhd(2,2), PCNhd(1,1)];

                            flag=1; %trigger the cell movement
                            k=1;    %count for number of cells moved so far in while loop
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12) % if the angle between current cell and cell being moved to is sufficiently shallow
                                    DirV(k)=DirV(k)-1;                          %remove one slab from current cell
                                    DirV(k+1)=DirV(k+1)+1;                      %add that cell to neighbor
                                    k=k+1;                                      %add one to cell moved count
                                    flag=1*(k<=dist);                           %repeat until k is equal to distance calculated earlier
            %next line for re-calculating chance of moving from new cell - if we wish to poll for plant cover at each step at some later time
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;                                     %if too steep an angle is encountered, movement stops 
                                end
                            end
                            Nhd(4,4)=DirV(1);                                   %replace neighborhood values with the direction vector values
                            Nhd(3,3)=DirV(2);
                            Nhd(2,2)=DirV(3);
                            Nhd(1,1)=DirV(4);

                            %NORTH - N
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (i>1)   %move with the wind / North
                            DirV=[Nhd(4,4), Nhd(3,4), Nhd(2,4), Nhd(1,4)];
                            PCNE=[PCNhd(4,4), PCNhd(3,4), PCNhd(2,4), PCNhd(1,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,4)=DirV(2);
                            Nhd(2,4)=DirV(3);
                            Nhd(1,4)=DirV(4);
                         
                            
                            %NORTH - NE
                        elseif Rand_dir>=((1+P)/2) && (i>1) && (j<size(Hstar,2))       %move to the starboard slab
                            DirV=[Nhd(4,4), Nhd(3,5), Nhd(2,6), Nhd(1,7)];
                            PCNE=[PCNhd(4,4), PCNhd(3,5), PCNhd(2,6), PCNhd(1,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,5)=DirV(2);
                            Nhd(2,6)=DirV(3);
                            Nhd(1,7)=DirV(4);
                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %NORTH - EAST BLOCK
                    
                    elseif WindDir==2 && rand<MoveChance      
                        
                        %MOVE N
                        if Rand_dir<((1-P)/2) %&& (i>1)
                            DirV=[Nhd(4,4), Nhd(3,4), Nhd(2,4), Nhd(1,4)];
                            PCNE=[PCNhd(4,4), PCNhd(3,4), PCNhd(2,4), PCNhd(1,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,4)=DirV(2);
                            Nhd(2,4)=DirV(3);
                            Nhd(1,4)=DirV(4);

%                             %MOVE NE
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2)% && (i>1) && (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(3,5), Nhd(2,6), Nhd(1,7)];
                            PCNE=[PCNhd(4,4), PCNhd(3,5), PCNhd(2,6), PCNhd(1,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,5)=DirV(2);
                            Nhd(2,6)=DirV(3);
                            Nhd(1,7)=DirV(4);

%                             %MOVE E
                        elseif Rand_dir>=((1+P)/2)% && (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(4,5), Nhd(4,6), Nhd(4,7)];
                            PCNE=[PCNhd(4,4), PCNhd(4,5), PCNhd(4,6), PCNhd(4,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,5)=DirV(2);
                            Nhd(4,6)=DirV(3);
                            Nhd(4,7)=DirV(4);
                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %EAST BLOCK
                    
                    elseif WindDir==3 && rand<MoveChance  
                        
                        %MOVES NE
                        if Rand_dir<((1-P)/2)% && (i>1) && (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(3,5), Nhd(2,6), Nhd(1,7)];
                            PCNE=[PCNhd(4,4), PCNhd(3,5), PCNhd(2,6), PCNhd(1,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,5)=DirV(2);
                            Nhd(2,6)=DirV(3);
                            Nhd(1,7)=DirV(4);

                            %MOVES E
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(4,5), Nhd(4,6), Nhd(4,7)];
                            PCNE=[PCNhd(4,4), PCNhd(4,5), PCNhd(4,6), PCNhd(4,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,5)=DirV(2);
                            Nhd(4,6)=DirV(3);
                            Nhd(4,7)=DirV(4);


                            %MOVES SE
                        elseif Rand_dir>=((1+P)/2)% && (i<size(Hstar,1)) && (j<size(Hstar,2)) 
                            DirV=[Nhd(4,4), Nhd(5,5), Nhd(6,6), Nhd(7,7)];
                            PCNE=[PCNhd(4,4), PCNhd(5,5), PCNhd(6,6), PCNhd(7,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,5)=DirV(2);
                            Nhd(6,6)=DirV(3);
                            Nhd(7,7)=DirV(4);
                            

                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %SOUTH-EAST BLOCK
                        
                    elseif WindDir==4 && rand<MoveChance     
                        
                        %MOVES E
                        if Rand_dir<((1-P)/2)% && (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(4,5), Nhd(4,6), Nhd(4,7)];
                            PCNE=[PCNhd(4,4), PCNhd(4,5), PCNhd(4,6), PCNhd(4,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,5)=DirV(2);
                            Nhd(4,6)=DirV(3);
                            Nhd(4,7)=DirV(4);


                            %MOVES SE
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2)% && (i<size(Hstar,1)) && (j<size(Hstar,2)) 
                            DirV=[Nhd(4,4), Nhd(5,5), Nhd(6,6), Nhd(7,7)];
                            PCNE=[PCNhd(4,4), PCNhd(5,5), PCNhd(6,6), PCNhd(7,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,5)=DirV(2);
                            Nhd(6,6)=DirV(3);
                            Nhd(7,7)=DirV(4);

                            %MOVES S
                        elseif Rand_dir>=((1+P)/2) && (i<size(Hstar,1))
                            DirV=[Nhd(4,4), Nhd(5,4), Nhd(6,4), Nhd(7,4)];
                            PCNE=[PCNhd(4,4), PCNhd(5,4), PCNhd(6,4), PCNhd(7,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,4)=DirV(2);
                            Nhd(6,4)=DirV(3);
                            Nhd(7,4)=DirV(4);

                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %SOUTH BLOCK
                        
                    elseif WindDir==5 && rand<MoveChance 
                        
                        %MOVES SE
                        if Rand_dir<((1-P)/2) %&& (i<size(Hstar,1)) && (j<size(Hstar,2))
                            DirV=[Nhd(4,4), Nhd(5,5), Nhd(6,6), Nhd(7,7)];
                            PCNE=[PCNhd(4,4), PCNhd(5,5), PCNhd(6,6), PCNhd(7,7)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,5)=DirV(2);
                            Nhd(6,6)=DirV(3);
                            Nhd(7,7)=DirV(4);

                            %MOVES S
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (i<size(Hstar,1))
                            DirV=[Nhd(4,4), Nhd(5,4), Nhd(6,4), Nhd(7,4)];
                            PCNE=[PCNhd(4,4), PCNhd(5,4), PCNhd(6,4), PCNhd(7,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,4)=DirV(2);
                            Nhd(6,4)=DirV(3);
                            Nhd(7,4)=DirV(4);

                            %MOVES SW
                        elseif Rand_dir>=((1+P)/2) %&& (i<size(Hstar,1)) %&& (j>1)
                            DirV=[Nhd(4,4), Nhd(5,3), Nhd(6,2), Nhd(7,1)];
                            PCNE=[PCNhd(4,4), PCNhd(5,3), PCNhd(6,2), PCNhd(7,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,3)=DirV(2);
                            Nhd(6,2)=DirV(3);
                            Nhd(7,1)=DirV(4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %SOUTH-WEST BLOCK
                        
                    elseif WindDir==6 && rand<MoveChance
                        
                        %MOVES S
                        if Rand_dir<((1-P)/2) %&& (i<size(Hstar,1))  
                            DirV=[Nhd(4,4), Nhd(5,4), Nhd(6,4), Nhd(7,4)];
                            PCNE=[PCNhd(4,4), PCNhd(5,4), PCNhd(6,4), PCNhd(7,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,4)=DirV(2);
                            Nhd(6,4)=DirV(3);
                            Nhd(7,4)=DirV(4);


                            %MOVES SW
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (i<size(Hstar,1)) && (j>1)  
                            DirV=[Nhd(4,4), Nhd(5,3), Nhd(6,2), Nhd(7,1)];
                            PCNE=[PCNhd(4,4), PCNhd(5,3), PCNhd(6,2), PCNhd(7,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,3)=DirV(2);
                            Nhd(6,2)=DirV(3);
                            Nhd(7,1)=DirV(4);


                            %MOVES W
                        elseif Rand_dir>=((1+P)/2) %&& (j>1) 
                            DirV=[Nhd(4,4), Nhd(4,3), Nhd(4,2), Nhd(4,1)];
                            PCNE=[PCNhd(4,4), PCNhd(4,3), PCNhd(4,2), PCNhd(4,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,3)=DirV(2);
                            Nhd(4,2)=DirV(3);
                            Nhd(4,1)=DirV(4);

                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %WEST BLOCK
                        
                    elseif WindDir==7 && rand<MoveChance   
                        
                        %MOVES SW
                        if Rand_dir<((1-P)/2) && (i<size(Hstar,1)) %&& (j>1)
                            DirV=[Nhd(4,4), Nhd(5,3), Nhd(6,2), Nhd(7,1)];
                            PCNE=[PCNhd(4,4), PCNhd(5,3), PCNhd(6,2), PCNhd(7,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(5,3)=DirV(2);
                            Nhd(6,2)=DirV(3);
                            Nhd(7,1)=DirV(4);

                            %MOVES W
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (j>1) 
                            DirV=[Nhd(4,4), Nhd(4,3), Nhd(4,2), Nhd(4,1)];
                            PCNE=[PCNhd(4,4), PCNhd(4,3), PCNhd(4,2), PCNhd(4,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,3)=DirV(2);
                            Nhd(4,2)=DirV(3);
                            Nhd(4,1)=DirV(4);

                            %MOVES NW
                        elseif Rand_dir>=((1+P)/2) %&& (i>1) && (j>1)
                            DirV=[Nhd(4,4), Nhd(3,3), Nhd(2,2), Nhd(1,1)];
                            PCNE=[PCNhd(4,4), PCNhd(3,3), PCNhd(2,2), PCNhd(1,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,3)=DirV(2);
                            Nhd(2,2)=DirV(3);
                            Nhd(1,1)=DirV(4);

                        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        %NORTH-WEST
                        
                    elseif WindDir==8 && rand<MoveChance      
                        
                        %MOVES W
                        if Rand_dir<((1-P)/2) %&& (j>1)
                            DirV=[Nhd(4,4), Nhd(4,3), Nhd(4,2), Nhd(4,1)];
                            PCNE=[PCNhd(4,4), PCNhd(4,3), PCNhd(4,2), PCNhd(4,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(4,3)=DirV(2);
                            Nhd(4,2)=DirV(3);
                            Nhd(4,1)=DirV(4);

                            %MOVES NW
                        elseif Rand_dir>=((1-P)/2) && Rand_dir<((1+P)/2) %&& (i>1) && (j>1)
                            DirV=[Nhd(4,4), Nhd(3,3), Nhd(2,2), Nhd(1,1)];
                            PCNE=[PCNhd(4,4), PCNhd(3,3), PCNhd(2,2), PCNhd(1,1)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,3)=DirV(2);
                            Nhd(2,2)=DirV(3);
                            Nhd(1,1)=DirV(4);
                            
                            %MOVES N
                        elseif Rand_dir>=((1+P)/2) %&& (i>1) 
                            DirV=[Nhd(4,4), Nhd(3,4), Nhd(2,4), Nhd(1,4)];
                            PCNE=[PCNhd(4,4), PCNhd(3,4), PCNhd(2,4), PCNhd(1,4)];
                            flag=1;
                            k=1;
                            while flag==1
                                if atan(((DirV(k+1)-DirV(k))*delta)/L)<(pi/12)
                                    DirV(k)=DirV(k)-1;
                                    DirV(k+1)=DirV(k+1)+1;
%                                     cellcnt=cellcnt+1;
%                                     MoveChance=(1-PCNE(k+1))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));
                                    k=k+1;
                                    flag=1*(k<=dist);
                                    MoveChance=(1-PCNE(k))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin))*MvF;
                                    if rand >= MoveChance
                                        flag=0;
                                    end
                                else
                                    flag=0;
                                end
                            end
                            Nhd(4,4)=DirV(1);
                            Nhd(3,4)=DirV(2);
                            Nhd(2,4)=DirV(3);
                            Nhd(1,4)=DirV(4);
%                             
%                             flag=1;
%                             cnt=0;
%                             k=1;
%                             while flag==1
%                                 if atan(((Nhd(4-k,4)-Nhd(4-k+1,4))*delta)/L)<(pi/12)
%                                     Nhd(4-k+1,4)=Nhd(4-k+1,4)-1;%4-k+1=4-(k-1)
%                                     Nhd(4-k,4)=Nhd(4-k,4)+1;
%                                     MoveChance=(1-PCNhd(4-k,4))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));%(calculating, but not using yet)
%                                     k=k+1;
%                                     cnt=cnt+1;
%                                     flag=1*(k<=dist);%*(Nhd(4,4)>=0)
%                                 else 
%                                     flag=0;
%                                 end
%                             end
                        end
                        end
                    end
            Hp(i-3:i+3,j-3:j+3,1)=Nhd;
            end
        else
        end
    end 
end

Hstar(:,:,1)=Hp(4:n1+3,4:n2+3);

end %end of AT function

%COMMENT: Old way of polling may have made it easier to change number of
%steps in each direction that a cell can possible move:

%it used this neighborhood function to declare a local region of 3 cells in
%each direction. If we used a variablelike D=#cells wished to move

%then...

% 
% function N1=NewNhd(R,C,H)
% if R<=(size(H,1)-1)
%     R=floor(R);
%     Rp=R+1;
% else
%     Rp=floor(R);
% end
% if R>=2
%     R=floor(R);
%     Rm=R-1;
% else
%     Rm=floor(R);
% end
% if C<=(size(H,2)-1)
%     C=floor(C);
%     Cp=C+1;
% else
%     Cp=floor(C);
% end
% if C>=2
%     C=floor(C);
%     Cm=C-1;
% else
%     Cm=floor(C);
% end
% N1=zeros(3,3);                %<------------CHANGED THIS TO zeros(D,D)******

% N1=[H(Rm,Cm),H(Rm,C),H(Rm,Cp);H(R,Cm),H(R,C),H(R,Cp);H(Rp,Cm),H(Rp,C),H(Rp,Cp)];
% 
% end

% *******************and used the old method of checking:

% *******************All of these 4's would need to be replaced with D+1

%                             flag=1;
%                             cnt=0;
%                             k=1;
%                             while flag==1
%                                 if atan(((Nhd(4-k,4)-Nhd(4-k+1,4))*delta)/L)<(pi/12)
%                                     Nhd(4-k+1,4)=Nhd(4-k+1,4)-1;%4-k+1=4-(k-1)
%                                     Nhd(4-k,4)=Nhd(4-k,4)+1;
%                                     MoveChance=(1-PCNhd(4-k,4))*((max(Windspeed-Windmin,0))/(StormThreshold-Windmin));%(calculating, but not using yet)
%                                     k=k+1;
%                                     cnt=cnt+1;
%                                     flag=1*(k<=dist);%*(Nhd(4,4)>=0)
%                                 else 
%                                     flag=0;
%                                 end
%                             end
