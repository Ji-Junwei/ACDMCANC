%% Distributed MANC depending on gradient tranmission
% 1 reference, K secondary source, K error sensors


classdef DMANC_CompensateSP
    properties
        Wc     % global control filter (K x (wlen+clen-1))
        Nabla  % local gradient (K x wlen)
        wlen   % length of control filter
        SecP   % secondary path estimates (K x K x slen)
        slen  % length of secondary path
        Xd    % Disturbance
        yc    % control signal
        Numnode  % number of nodes, K
        C     % compensate filter (M x M x clen)
        clen  % length of compensate filter
    end


    methods
        % initial
        function obj = DMANC_CompensateSP(wLen,SecondaryPath,sLen,node_num,N,Dis,cLen)
            obj.wlen    = wLen;
            obj.Wc      = zeros(node_num,(wLen));
            obj.Nabla   = zeros(node_num,(wLen+cLen-1));
            obj.SecP    = SecondaryPath;
            obj.slen    = sLen;
            obj.Numnode = node_num;
            obj.yc      = zeros(node_num,N);
            obj.Xd = Dis;
            obj.C       = zeros(node_num,node_num,cLen);
            obj.clen    = cLen;
        end

       % obtain compensation filter
        function [err,obj] = CompensateSP(obj,muc)
            T   = 200000;           % duration
            wgn = randn(1,T);       % generate white noise
            err = zeros(obj.Numnode,obj.Numnode,T); % error signal

            for m = 1: obj.Numnode
                for k = 1: obj.Numnode
                    if m == k
                        continue;
                    else
                        wgnd = filter(reshape(obj.SecP(m,k,:),[1,obj.slen]),1,wgn);                     % wgn pass cross secondary path
                        FwgnLMS = dsp.FilteredXLMSFilter(obj.clen,"StepSize",muc,...
                        "SecondaryPathCoefficients",reshape(obj.SecP(m,m,:),[1,obj.slen]),...
                        "SecondaryPathEstimate",reshape(obj.SecP(m,m,:),[1,obj.slen]));
                        [~,e] = FwgnLMS(wgn,wgnd);
                        err(m,k,:) = e;
                        CF = -flip(FwgnLMS.Coefficients);
                        obj.C(m,k,:) = reshape(CF,[1,1,obj.clen]);
                    end
                end
            end

        end


        % MGDFxLMS; 166
        function [e,obj] = DMANC_gradient_166(obj,xin,muw)
            N = length(xin);       %duration
            e = zeros(obj.Numnode,N);       %error signal K x duration
            xc = zeros(1,obj.wlen);      % reference vector


            ys = zeros(obj.Numnode,obj.slen);               % y buffer for filter secondary path
            xs = zeros(1,obj.slen);      % filter secondary path
            xf = zeros(obj.Numnode,(obj.wlen+obj.clen-1));     % filtered reference

            for i =1:N
                xc = [xin(i) xc(1:(end-1))];   % update reference vector
                

                % generate control signal
                obj.yc(1,i) = sum(obj.Wc(1,:).*xc);
                obj.yc(2,i) = sum(obj.Wc(2,:).*xc);
                obj.yc(3,i) = sum(obj.Wc(3,:).*xc);
                obj.yc(4,i) = sum(obj.Wc(4,:).*xc);
                obj.yc(5,i) = sum(obj.Wc(5,:).*xc);
                obj.yc(6,i) = sum(obj.Wc(6,:).*xc);
                
                % error sensors received
                ys(1,:) = [obj.yc(1,i) ys(1,1:(obj.slen-1))];                            % y1 buffer update
                ys(2,:) = [obj.yc(2,i) ys(2,1:(obj.slen-1))];                            % y2 buffer update
                ys(3,:) = [obj.yc(3,i) ys(3,1:(obj.slen-1))];                            % y3 buffer update
                ys(4,:) = [obj.yc(4,i) ys(4,1:(obj.slen-1))];                            % y4 buffer update
                ys(5,:) = [obj.yc(5,i) ys(5,1:(obj.slen-1))];                            % y5 buffer update
                ys(6,:) = [obj.yc(6,i) ys(6,1:(obj.slen-1))];                            % y6 buffer update

                % error 1
                y1 = sum(reshape(obj.SecP(1,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(1,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(1,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(1,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(1,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(1,6,:),[1,obj.slen]).*ys(6,:));
                e(1,i) = obj.Xd(1,i)-y1;
                % error 2
                y2 = sum(reshape(obj.SecP(2,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(2,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(2,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(2,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(2,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(2,6,:),[1,obj.slen]).*ys(6,:));
                e(2,i) = obj.Xd(2,i)-y2;
                % error 3
                y3 = sum(reshape(obj.SecP(3,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(3,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(3,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(3,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(3,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(3,6,:),[1,obj.slen]).*ys(6,:));
                e(3,i) = obj.Xd(3,i)-y3;
                % error 4
                y4 = sum(reshape(obj.SecP(4,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(4,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(4,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(4,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(4,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(4,6,:),[1,obj.slen]).*ys(6,:));
                e(4,i) = obj.Xd(4,i)-y4;
                % error 5
                y5 = sum(reshape(obj.SecP(5,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(5,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(5,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(5,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(5,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(5,6,:),[1,obj.slen]).*ys(6,:));
                e(5,i) = obj.Xd(5,i)-y5;
                % error 6
                y6 = sum(reshape(obj.SecP(6,1,:),[1,obj.slen]).*ys(1,:)) + sum(reshape(obj.SecP(6,2,:),[1,obj.slen]).*ys(2,:)) + ...
                     sum(reshape(obj.SecP(6,3,:),[1,obj.slen]).*ys(3,:)) + sum(reshape(obj.SecP(6,4,:),[1,obj.slen]).*ys(4,:)) + ...
                     sum(reshape(obj.SecP(6,5,:),[1,obj.slen]).*ys(5,:)) + sum(reshape(obj.SecP(6,6,:),[1,obj.slen]).*ys(6,:));
                e(6,i) = obj.Xd(6,i)-y6;

                % gradient
                xs = [xin(i) xs(1:(end-1))];

                % controller 1
                fx1 = sum(xs.*reshape(obj.SecP(1,1,:),[1,obj.slen]));
                xf(1,:) = [fx1 xf(1,1:(end-1))];
                obj.Nabla(1,:) = muw*xf(1,:)*e(1,i);

                % controller 2
                fx2 = sum(xs.*reshape(obj.SecP(2,2,:),[1,obj.slen]));
                xf(2,:) = [fx2 xf(2,1:(end-1))];
                obj.Nabla(2,:) = muw*xf(2,:)*e(2,i);

                % controller 3
                fx3 = sum(xs.*reshape(obj.SecP(3,3,:),[1,obj.slen]));
                xf(3,:) = [fx3 xf(3,1:(end-1))];
                obj.Nabla(3,:) = muw*xf(3,:)*e(3,i);

                % controller 4
                fx4 = sum(xs.*reshape(obj.SecP(4,4,:),[1,obj.slen]));
                xf(4,:) = [fx4 xf(4,1:(end-1))];
                obj.Nabla(4,:) = muw*xf(4,:)*e(4,i);

                % controller 5
                fx5 = sum(xs.*reshape(obj.SecP(5,5,:),[1,obj.slen]));
                xf(5,:) = [fx5 xf(5,1:(end-1))];
                obj.Nabla(5,:) = muw*xf(5,:)*e(5,i);

                % controller 6
                fx6 = sum(xs.*reshape(obj.SecP(6,6,:),[1,obj.slen]));
                xf(6,:) = [fx6 xf(6,1:(end-1))];
                obj.Nabla(6,:) = muw*xf(6,:)*e(6,i);


                % generate global control filter
                % controller 1
                  % filtered compensation filter
                  a1 = filter(reshape(obj.C(2,1,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a2 = filter(reshape(obj.C(3,1,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,1,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,1,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,1,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(1,:) = obj.Wc(1,:) + obj.Nabla(1,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 2
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,2,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(3,2,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a3 = filter(reshape(obj.C(4,2,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,2,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,2,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(2,:) = obj.Wc(2,:) + obj.Nabla(2,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);
                
                % controller 3
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,3,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,3,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(4,3,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a4 = filter(reshape(obj.C(5,3,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,3,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(3,:) = obj.Wc(3,:) + obj.Nabla(3,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);
                  
                
                % controller 4
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,4,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,4,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,4,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(5,4,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  a5 = filter(reshape(obj.C(6,4,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(4,:) = obj.Wc(4,:) + obj.Nabla(4,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 5
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,5,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,5,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,5,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,5,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a5 = filter(reshape(obj.C(6,5,:),[1,obj.clen]),1,obj.Nabla(6,:));
                  obj.Wc(5,:) = obj.Wc(5,:) + obj.Nabla(5,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

                % controller 6
                % filtered compensation filter
                  a1 = filter(reshape(obj.C(1,6,:),[1,obj.clen]),1,obj.Nabla(1,:));
                  a2 = filter(reshape(obj.C(2,6,:),[1,obj.clen]),1,obj.Nabla(2,:));
                  a3 = filter(reshape(obj.C(3,6,:),[1,obj.clen]),1,obj.Nabla(3,:));
                  a4 = filter(reshape(obj.C(4,6,:),[1,obj.clen]),1,obj.Nabla(4,:));
                  a5 = filter(reshape(obj.C(5,6,:),[1,obj.clen]),1,obj.Nabla(5,:));
                  obj.Wc(6,:) = obj.Wc(6,:) + obj.Nabla(6,1:obj.wlen) + a1((end-obj.wlen+1):end) +...
                                              a2((end-obj.wlen+1):end) + a3((end-obj.wlen+1):end) + ...
                                              a4((end-obj.wlen+1):end) + a5((end-obj.wlen+1):end);

            end


        end

       

    end
end