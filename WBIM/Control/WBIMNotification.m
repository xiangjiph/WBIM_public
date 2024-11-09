classdef WBIMNotification < handle
    
    properties(Constant, Access=private)
        SENDER_EMAIL = 'hello_reviewer@gmail.com';
        SENDER_PASSCODE = 'hello_reviewer';
        RECEIVER_EMAIL = {'hello_reviewer@ucsd.edu', 'hello_reviewer@ucsd.edu'}; 
    end
    
    methods
        function obj = WBIMNotification()
            obj.initialization();
        end
        
        function send_email(obj, tilte, message)
            try
                message_head = sprintf('%s', datestr(now, 'yyyy-mm-dd hh:MM:ss.FFF'));
                if nargin < 3
                    message = message_head;
                else
                    message = sprintf('%s\t%s', message_head, message);
                end
                for i = 1:length(obj.RECEIVER_EMAIL)
                    sendmail(obj.RECEIVER_EMAIL{i}, tilte, message);
                end
            catch ME
                warning("testFailed to send the email. Error message:\n%s", ...\
                    getReport(ME, 'extended', 'hyperlinks', 'off'));
            end
        end
    end
    
    methods(Access=private)
        
        function initialization(obj, overwriteQ)
            if nargin < 2
                overwriteQ = false;
            end
            persistent initialized_Q;
            if isempty(initialized_Q) || ~initialized_Q || overwriteQ
                %set up SMTP service for Gmail
                setpref('Internet','E_mail',obj.SENDER_EMAIL);
                setpref('Internet','SMTP_Server','smtp.gmail.com');
                setpref('Internet','SMTP_Username', obj.SENDER_EMAIL);
                setpref('Internet','SMTP_Password', obj.SENDER_PASSCODE);
                % Gmail server.
                props = java.lang.System.getProperties;
                props.setProperty('mail.smtp.auth','true');
                props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
                props.setProperty('mail.smtp.socketFactory.port','465');
                initialized_Q = true;
                fprintf('Finish initializing WBIMNotification\n');
            end
        end
        
    end
end