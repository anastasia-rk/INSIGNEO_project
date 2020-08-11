%%

%% Field decomposition
% figure; set(gcf,'color','w');
% imshow(A); hold on;
% side1 = knots(1,2) - knots(1,1);
% side2 = knots(2,2) - knots(2,1);
%  colors = ['g', 'y', 'w', 'r', 'c','m','g', 'y', 'w', 'r', 'c','m','g', 'y', 'w', 'r', 'c','m','g', 'y', 'w', 'r', 'c','m'];
% for i=1:length(Theta)
%     a = knots(1,i*2-1);
%     b = knots(2,i*2-1);
%     p=rectangle('Position',[a,b,side1,side2],'Curvature',0.1,'EdgeColor',colors(i)); 
%     hold on;
%     if i==i
%     txt = (['$\beta_{' num2str(i) '}$']);
%     xx = a + side1/2 - 60;
%     yy = b + side2/2 - 30;
%     text(xx,yy, 2,txt,'Color',colors(i),'FontSize',24)
%     end
% end
% ChangeInterpreter(gcf,'Latex');
% tightfig

%% Plot b-spline grid in 3d
% x0 = 20;
% y0 = 20;
% width=700;
% height=500;
% figure; set(gcf,'color','w'); 
% set(gcf,'units','points','position',[x0,y0,width,height]);
% image(A); hold on;
% view(3)
% z_ticklabels{1} = '';
% k = 1;
% z_ticklabels{k+1} = num2str(k-1);
% for j=1:4 %    [4 8 12 16 20]
%     a = knots(1,j*2-1);
%     b = knots(2,j*2-1);
%     Theta_temp = Theta;
%     Theta_temp(j) = k;
%     k = k + 1;
%     z_ticklabels{k+1} = num2str(k-1);
%     gridlim = [a b a+side1 b+side2];
%     done = plot_surface(Theta_temp,Z,knots,grid_limits,basis_type);
%     hold on;
%     txt = (['$\beta_{' num2str(j) '}$']);
%     xx = a + side1/2 - 30;
%     yy = b + side2/2 - 30;
%     text(xx,yy,k-0.5,txt,'Color','k','FontSize',24)
% end
% % view(3)
% az = -30;
% el = 40;
% view(az, el);
% zlim([0 k]);
% z_ticks = [1:1:k];
% set(gca,'ZTick',[0 z_ticks])
% set(gca,'XTick',[])
% set(gca,'YTick',[])
% set(gca,'ZTickLabel',z_ticklabels)
% zlabel('Grid of basis functions')
% ylabel('Brightfield')
% ChangeInterpreter(gcf,'Latex');
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
% %%
% tracks = 64;
% 
% figure; set(gcf,'color','w');
% imshow(Cnt); hold on;
% for j = tracks
%    plot(Y{j}(:,1)-padW,       Y{j}(:,2)-padH,'-g','LineWidth',2); hold on; 
% end
% cleanfigure;
%  matlab2tikz('track_fish.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '6cm', 'width','6cm');
% 
% %%
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(          Y{j}(:,1)./cc,       Y{j}(:,2)./cc,'-k','LineWidth',1.5); hold on; 
%    plot(X_cond{j,1}(1,:)./cc,X_cond{j,1}(2,:)./cc,'--m','LineWidth',1.5); hold on;
%    plot(X_cond{j,2}(1,:)./cc,X_cond{j,2}(2,:)./cc,'--c','LineWidth',1.5); hold on;
%    plot(      X_out{j}(1,:)./cc,      X_out{j}(2,:)./cc,':g','LineWidth',1.5); hold on;
%    plot(      X_out{j}(1,1)./cc,      X_out{j}(2,1)./cc,'og','LineWidth',1.5); hold on;
% end
% xlabel('X, px'); ylabel('Y, px');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
%  cleanfigure;
%  matlab2tikz('urts_track.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '6cm', 'width','6cm');
% 
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond_f{j,1}(1,end:2)./cc,'--m','LineWidth',1.5); hold on;
%    plot(X_cond_f{j,2}(1,end:2)./cc,'.c','LineWidth',1.5); hold on;
%    plot(X_filt{j}(1,2:end)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_filt{j}(1,2:end)./cc,'.g','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$s_x$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('ukf_sx.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond{j,1}(1,:)./cc,'--r','LineWidth',1.5); hold on;
%    plot(X_cond{j,2}(1,:)./cc,'.b','LineWidth',1.5); hold on;
%    plot(X_out{j}(1,:)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_out{j}(1,1)./cc,'og','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$s_x$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('urts_sx.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
%     
%      
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond_f{j,1}(2,end:2)./cc,'--m','LineWidth',1.5); hold on;
%    plot(X_cond_f{j,2}(2,end:2)./cc,'.c','LineWidth',1.5); hold on;
%    plot(X_filt{j}(2,2:end)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_filt{j}(2,2:end)./cc,'.g','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$s_y$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('ukf_sy.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond{j,1}(2,:)./cc,'--r','LineWidth',1.5); hold on;
%    plot(X_cond{j,2}(2,:)./cc,'.b','LineWidth',1.5); hold on;
%    plot(X_out{j}(2,:)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_out{j}(2,1)./cc,'og','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$s_y$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('urts_sy.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
%      
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond_f{j,1}(3,end:2)./cc,'--m','LineWidth',1.5); hold on;
%    plot(X_cond_f{j,2}(3,end:2)./cc,'.c','LineWidth',1.5); hold on;
%    plot(X_filt{j}(3,2:end)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_filt{j}(3,2:end)./cc,'.g','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$v_x$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('ukf_vx.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond{j,1}(3,:)./cc,'--r','LineWidth',1.5); hold on;
%    plot(X_cond{j,2}(3,:)./cc,'.b','LineWidth',1.5); hold on;
%    plot(X_out{j}(3,:)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_out{j}(3,1)./cc,'og','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$v_x$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('urts_vx.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
%     
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond_f{j,1}(4,end:2)./cc,'--m','LineWidth',1.5); hold on;
%    plot(X_cond_f{j,2}(4,end:2)./cc,'.c','LineWidth',1.5); hold on;
%    plot(X_filt{j}(4,2:end)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_filt{j}(4,2:end)./cc,'.g','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$v_y$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('ukf_vy.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
% figure; set(gcf,'color','w','OuterPosition', [100, 100, 500, 300]);
% for j = tracks
%    plot(X_cond{j,1}(4,:)./cc,'--r','LineWidth',1.5); hold on;
%    plot(X_cond{j,2}(4,:)./cc,'.b','LineWidth',1.5); hold on;
%    plot(X_out{j}(4,:)./cc,'-g','LineWidth',1.5); hold on;
%    plot(X_out{j}(4,1)./cc,'og','LineWidth',1.5); hold on;
% end
% xlim([0,length(X_out{j})]);
% ylabel('$v_y$, $\mu$m'); xlabel('t, min');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% ChangeInterpreter(gcf,'Latex');
% cleanfigure;
% matlab2tikz('urts_vy.tikz', 'showInfo', false,'parseStrings',false, ...
%          'standalone', false,'height', '3cm', 'width','4cm');
    


%% Convergence plots

% fig8 
% parameter convergence
% figure; set(gcf,'color','w');
% plot(DT_plot,'ko','Linewidth',1);
% ylabel('$\Delta\theta$', 'interpreter', 'latex');
% xlabel('\textrm{Iteration}');
% ChangeInterpreter(gcf,'Latex');

% fig8 
% Q convergence
% figure; set(gcf,'color','w');
% plot(LL_plot,'ko','Linewidth',1);
% ylabel('$\mathcal{Q}$', 'interpreter', 'latex');
% xlabel('\textrm{Iteration}');
% ChangeInterpreter(gcf,'Latex');
% 
% figure; set(gcf,'color','w');
% plot(LL_plot,'k','Linewidth',2);
% ylabel('$\mathcal{Q}$', 'interpreter', 'latex');
% xlabel('\textrm{Iteration}');
% ChangeInterpreter(gcf,'Latex');



% figure; set(gcf,'color','w');
% subplot(5,1,1);
% Track1 = 10;
% stairs(Mode{Track1},'k','Linewidth',1);
% ylim([1;3]); xlim([1,inf]); xlabel('time, sec'); ylabel('mode');
% subplot(5,1,2);
% Track2 = 36;
% stairs(Mode{Track2},'k','Linewidth',1);
% ylim([1;3]); xlim([1,inf]); xlabel('time, sec'); ylabel('mode');
% subplot(5,1,3);
% Track3 = 65;
% stairs(Mode{Track3},'k','Linewidth',1);
% ylim([1;3]); xlim([1,inf]); xlabel('time, sec'); ylabel('mode');
% subplot(5,1,4);
% Track4 = 99;
% stairs(Mode{Track4},'k','Linewidth',1);
% ylim([1;3]); xlim([1,inf]); xlabel('time, sec'); ylabel('mode');
% subplot(5,1,5);
% Track5 = 168;
% stairs(Mode{Track5},'k','Linewidth',1);
% ylim([1;3]); xlim([1,inf]); xlabel('time, sec'); ylabel('mode');
% ChangeInterpreter(gcf,'Latex');

% %% Track
% figure; set(gcf,'color','w');
% j = 36;
% imshow(A); hold on;
% plot(X_rts{j}(1:end-1,1),X_rts{j}(1:end-1,2),'-w','LineWidth',1); hold on;
% if ind_resample{j}~=0
% for jj = ind_resample{j}
%     plot(X_resample{jj}(1:end-1,1),X_resample{jj}(1:end-1,2),'-g','LineWidth',1); hold on;
% end
% end
% if ind_dead{j}~=0
%   for jj = ind_dead{j}
%    plot(X_dead{jj}(:,1),X_dead{jj}(:,2),'-r','LineWidth',1); hold on; 
%   end
% end

% %% Distance from wound
% figure;set(gcf,'color','w');
% j = 11;
% plot(990 - X_rts{j}(1:end-1,1),'-k','LineWidth',1); hold on;
% if ind_resample{j}~=0
% for jj = ind_resample{j}
%     plot(990 - X_resample{jj}(1:end-1,1),'-g','LineWidth',1); hold on;
% end
% end
% if ind_dead{j}~=0
%   for jj = ind_dead{j}
%    plot(990 - X_dead{jj}(:,1),'-r','LineWidth',1); hold on; 
%   end
% end

%% Save figures to latex format
% h = figure(2); 
% h.Units = 'centimeters'; % set figure position to cm
% set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% ChangeInterpreter(gcf,'Latex');
% Plot2LaTeX( h, 'cell_modes' ) 
% 
% h = figure(2); 
% h.Units = 'centimeters'; % set figure position to cm
% set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% ChangeInterpreter(gcf,'Latex');
% Plot2LaTeX( h, 'cell_modes' ) 
% 
% h = figure(3); 
% h.Units = 'centimeters'; % set figure position to cm
% %h.Position(2) = [h.Position(2)-9]; % set figure position bevore resize
% % h.Position([3:4]) = [20,12]; % resize figure
% set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% ChangeInterpreter(gcf,'Latex');
% Plot2LaTeX( h, 'gradient' )
% 
% h = figure(4); 
% h.Units = 'centimeters'; % set figure position to cm
% %h.Position(2) = [h.Position(2)-9]; % set figure position bevore resize
% % h.Position([3:4]) = [20,12]; % resize figure
% set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% ChangeInterpreter(gcf,'Latex');
% Plot2LaTeX( h, 'smoother' )
% 
% h = figure(2); 
% h.Units = 'centimeters'; % set figure position to cm
% set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% ChangeInterpreter(gcf,'Latex');
% imagename = ['tracks_sigmas03' num2str(sig2_Q) '01'];
% imagename = ['tracks_sigmas010301'];
% Plot2LaTeX( h, imagename)
% % 
% % 
% 
% h = figure(4); 
% h.Units = 'centimeters'; % set figure position to cm
% % set(findall(h,'-property','FontSize'),'FontSize',10); % change font size
% % ChangeInterpreter(gcf,'Latex');
% imagename = ['trynow'];
% Plot2LaTeX( h, imagename)

% Tr = [1, 44, 62, 17, 111, 16, 38, 19, 64, 119 124, 88];
% 
% figure; set(gcf,'color','w');
% for j = 1:length(Tr)
%    subplot(4,3,j,'align')
%    plot(Y{Tr(j)}(:,1),Y{Tr(j)}(:,2),'-k','LineWidth',2); hold on; 
%    plot(Y{Tr(j)}(1,1),Y{Tr(j)}(1,2),'*k','LineWidth',1); hold on; grid on;
% %    xlabel('\textrm{X,px}', 'interpreter', 'latex');
% %    ylabel('\textrm{Y,px}', 'interpreter', 'latex');
%    title(['Track ' num2str(Tr(j))]);
% end
% hold off;
% ChangeInterpreter(gcf,'Latex');
% %% Twists and turns
% Tr = [19, 62, 44, 30, 124, 53, 16, 38, 1];
% 
% figure; set(gcf,'color','w');
% for j = 1:length(Tr)
%    subplot(3,3,j,'align')
% %    plot(Y{Tr(j)}(:,1),Y{Tr(j)}(:,2),'-k','LineWidth',1); hold on; 
% %    plot(Y{Tr(j)}(1,1),Y{Tr(j)}(1,2),'*k','LineWidth',1); hold on; grid on;
%    plot(X_rts{Tr(j)}(:,1),X_rts{Tr(j)}(:,2),'LineWidth',1); hold on; 
%    plot(X_rts{Tr(j)}(1,1),X_rts{Tr(j)}(1,2),'o','LineWidth',1); hold on;
% %    xlabel('\textrm{X,px}', 'interpreter', 'latex');
% %    ylabel('\textrm{Y,px}', 'interpreter', 'latex');
%    m1 = mod(j,3);
%    m2 = floor(j/3) + 1;
%    if j == 9
%        m2 = 3;
%    end
%    title(['Track ' num2str(Tr(j)) ': $ M^{' num2str(m1) '} \rightarrow M^{' num2str(m2) '}$' ]);
% end
% hold off;
% set(findall(gcf,'-property','FontSize'),'FontSize',18);
% ChangeInterpreter(gcf,'Latex');
% 
% tightfig

% %% Selected tracks
% gg = [0.5 0.5 0.5];
% figure; set(gcf,'color','w');
%  imshow(A); hold on;
%  alpha(0.6);
% for j = 1:nTracks
%    plot(Y{j}(:,1),Y{j}(:,2),'Color',gg,'LineWidth',1); hold on; 
% end
% for j = 1:length(Tr)
%    i = Tr(j);
%    plot(Y{i}(:,1),Y{i}(:,2),'LineWidth',2); hold on; 
%    tx = num2str(i);
%    text(Y{i}(1,1),Y{i}(1,2),2, tx, 'interpreter', 'latex','Color','k','FontSize',18);
% end
% surf(Yy_grid,Xx_grid,-AA,'FaceColor',white,'EdgeColor',white);
% view(2)
% line([x_max-75-padW,x_max-75-padW],[210,770],'Color','k','LineWidth',2);
% txt1 = (' Wound at 600 $\mu$m');
% text(x_max - 300,210, 2,txt1,'Color','k','FontSize',14)
% hold on;
% % line([50,221],[575,575],'Color','k','LineWidth',5);
% line([50,221],[800,800],'Color','k','LineWidth',5);
% txt = ('100 $\mu$m');
% text(80,780, 2,txt,'Color','k','FontSize',14)
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
% ChangeInterpreter(gcf,'Latex');
% 
% tightfig

%% Making tikz
% cleanfigure;
% matlab2tikz('example_tracks.tikz', 'showInfo', false, ...
%         'parseStrings',false,'standalone', false, ...
%         'height', '\figureheight', 'width','\figurewidth');