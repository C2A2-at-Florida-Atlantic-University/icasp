classdef utils
    methods (Static)
        % Calculates ground truth DoA angles, based on client coordinates (assumes that RX is at 0;0;0)
        function [azimuth, elevation] = get_ground_truth(client_coord)
            azimuth = round(rad2deg(atan(client_coord(1) / client_coord(2))), 4) * -1;
            elevation = round(rad2deg(atan(sqrt(client_coord(1)^2 + client_coord(2)^2) / client_coord(3))), 4);
        end

        function [frame_corr] = run_doppler_calibration(frame, frame_len, lts_start, fig_title, lts_len, show_figures)
            % Define constants
            LTS1_START_IDX = lts_start;
            LTS2_START_IDX = LTS1_START_IDX + lts_len;
            LTS3_START_IDX = LTS2_START_IDX + lts_len;
        
            % Extract LTS1 & LTS2 IQ samples
            ant1_lts1 = frame(LTS1_START_IDX : LTS2_START_IDX - 1);
            ant1_lts2 = frame(LTS2_START_IDX : LTS3_START_IDX - 1);
            
            % Calculate Doppler shift between LTS1 & LTS2
            w12 = angle(ant1_lts1'*ant1_lts2)/64;
            
            % Apply correction to the whole frame
            frame_corr = frame.*exp(-1j * w12 * [0:frame_len-1]');
        
            % Prepare corrected LTS1 & LTS2 for calibration
            ant1_lts1_corr = frame_corr(LTS1_START_IDX : LTS2_START_IDX - 1);
            ant1_lts2_corr = frame_corr(LTS2_START_IDX : LTS3_START_IDX - 1);
            
            % Plot
            if show_figures
                figure
                
                subplot(2,1,1)
                plot(1:lts_len, ant1_lts1)
                hold on
                plot(1:lts_len, ant1_lts2)
                hold off
                title("Raw LTS1 & LTS2 sequences.")
            
                subplot(2,1,2)
                plot(1:lts_len, ant1_lts1_corr)
                hold on
                plot(1:lts_len, ant1_lts2_corr)
                hold off
                title("Corrected LTS1 & LTS2 sequences.")
            
                sgtitle(fig_title)
            end
        end

        function [ant1_final, ant2_final, ant3_final, ant4_final, w_ant12, w_ant13, w_ant14] = run_interant_calibration4(lts_indexes, ant1, ant2, ant3, ant4, LTS_LEN, LTS_COUNT)
            ant1_sub = ant1(lts_indexes(1) : lts_indexes(1) + LTS_LEN * LTS_COUNT - 1);
            ant2_sub = ant2(lts_indexes(2) : lts_indexes(2) + LTS_LEN * LTS_COUNT - 1);
            ant3_sub = ant3(lts_indexes(3) : lts_indexes(3) + LTS_LEN * LTS_COUNT - 1);
            ant4_sub = ant4(lts_indexes(4) : lts_indexes(4) + LTS_LEN * LTS_COUNT - 1);
        
            % Calculate phase offsets between these frames
            w_ant12 = angle(ant1_sub'*ant2_sub);
            w_ant13 = angle(ant1_sub'*ant3_sub);
            w_ant14 = angle(ant1_sub'*ant4_sub);
            
            % Apply corrections
            ant1_final = ant1;
            ant2_final = ant2.*exp(-1j * (w_ant12));
            ant3_final = ant3.*exp(-1j * (w_ant13));
            ant4_final = ant4.*exp(-1j * (w_ant14));
        end

        function [ant1_final, ant2_final, ant3_final, ant4_final, ant5_final, ant6_final, w_ant12, w_ant13, w_ant14, w_ant15, w_ant16] = run_interant_calibration6(lts_indexes, ant1, ant2, ant3, ant4, ant5, ant6, LTS_LEN, LTS_COUNT)
            ant1_sub = ant1(lts_indexes(1) : lts_indexes(1) + LTS_LEN * LTS_COUNT - 1);
            ant2_sub = ant2(lts_indexes(2) : lts_indexes(2) + LTS_LEN * LTS_COUNT - 1);
            ant3_sub = ant3(lts_indexes(3) : lts_indexes(3) + LTS_LEN * LTS_COUNT - 1);
            ant4_sub = ant4(lts_indexes(4) : lts_indexes(4) + LTS_LEN * LTS_COUNT - 1);
            ant5_sub = ant5(lts_indexes(5) : lts_indexes(5) + LTS_LEN * LTS_COUNT - 1);
            ant6_sub = ant6(lts_indexes(6) : lts_indexes(6) + LTS_LEN * LTS_COUNT - 1);
        
            % Calculate phase offsets between these frames
            w_ant12 = angle(ant1_sub'*ant2_sub);
            w_ant13 = angle(ant1_sub'*ant3_sub);
            w_ant14 = angle(ant1_sub'*ant4_sub);
            w_ant15 = angle(ant1_sub'*ant5_sub);
            w_ant16 = angle(ant1_sub'*ant6_sub);
            
            % Apply corrections
            ant1_final = ant1;
            ant2_final = ant2.*exp(-1j * (w_ant12));
            ant3_final = ant3.*exp(-1j * (w_ant13));
            ant4_final = ant4.*exp(-1j * (w_ant14));
            ant5_final = ant5.*exp(-1j * (w_ant15));
            ant6_final = ant6.*exp(-1j * (w_ant16));
        end

        function [iq_final, w_ant] = run_interant_calibrationN_rows(lts_positions, iq_grid, LTS_LEN, LTS_COUNT, ANT_ROWS)
            % Prepare 2D array to contain inter-antenna offsets (w.r.t. ant #1) 
            w_ant = zeros(size(iq_grid, 1), size(iq_grid, 2));

            % Prepare 3D array to contain calibrated IQ samples
            iq_final = zeros(size(iq_grid));

            % Run separate calibration for each row
            for i = 1:size(iq_grid, 1)
                % Retrieve LTS indexes for a row
                lts_indexes = lts_positions(ANT_ROWS(i, :)) + 1;

                % Extract antenna IQ samples
                ant1 = reshape(iq_grid(i, 1, :), [2048 1]);
                ant2 = reshape(iq_grid(i, 2, :), [2048 1]);
                ant3 = reshape(iq_grid(i, 3, :), [2048 1]);
                ant4 = reshape(iq_grid(i, 4, :), [2048 1]);

                % Perform calibration
                [ant1_final, ant2_final, ant3_final, ant4_final, w_ant12, w_ant13, w_ant14] = calibration.run_interant_calibration4(lts_indexes, ant1, ant2, ant3, ant4, LTS_LEN, LTS_COUNT);

                w_ant(i, 1) = 1;
                w_ant(i, 2) = w_ant12;
                w_ant(i, 3) = w_ant13;
                w_ant(i, 4) = w_ant14;

                iq_final(i, 1, :) = ant1_final;
                iq_final(i, 2, :) = ant2_final;
                iq_final(i, 3, :) = ant3_final;
                iq_final(i, 4, :) = ant4_final;
            end
        end

        function [iq_final, w_ant] = run_interant_calibrationN(lts_positions, iq_grid, LTS_LEN, LTS_COUNT, ANT_ROWS)
            % Prepare 2D array to contain inter-antenna offsets (w.r.t. ant #1) 
            w_ant = zeros(size(iq_grid, 1), size(iq_grid, 2));

            % Prepare 3D array to contain calibrated IQ samples
            iq_final = zeros(size(iq_grid));

            % Prepare 1st antenna samples
            ant1_lts_index = lts_positions(ANT_ROWS(1, 1)) + 1;
            ant1_sub = iq_grid(1, 1, ant1_lts_index : ant1_lts_index + LTS_LEN * LTS_COUNT - 1);

            iq_final(1, 1, :) = iq_grid(1, 1, :);

            for i = 1:size(w_ant, 1)
                for j = 1:size(w_ant, 2)
                    antk_lts_index = lts_positions(ANT_ROWS(i, j)) + 1;
                    antk_sub = iq_grid(i, j, antk_lts_index : antk_lts_index + LTS_LEN * LTS_COUNT - 1);

                    % Calculate phase offsets between antennas 1 and K
                    w_ant(i, j) = angle(ant1_sub(:)'*antk_sub(:));

                    % Apply corrections
                    iq_final(i, j, :) = iq_grid(i, j, :).*exp(-1j * (w_ant(i, j)));
                end
            end
        end
    end
end