#!/user/bin/env python3

import sys
import os
import argparse
import logging
import smtplib
from email.message import EmailMessage

from astropy.table import Table

'''
To check if a SBID was deposited/released
'''

# log_stream = StringIO()
log = logging.getLogger(__name__)

__author__ = "Yuanming Wang <yuanmingwang@swin.edu.au>"


def _main():
    parser = argparse.ArgumentParser(prog='VOevent', description='VO Event trigger')
    parser.add_argument('sbid', type=int, help='SBID in format of number')
    parser.add_argument('-p', '--path', type=str, default='.', 
                        help='path where SBxxxx were stored')
    parser.add_argument('-e', '--email', default=None, nargs='+',
                        help='send email notification to a list of email addresses')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='make it verbose')
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(
            # stream=log_stream, 
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.DEBUG,
            datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(
            # stream=log_stream, 
            format='%(asctime)s.%(msecs)03d %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')

    log.debug(args)

    fname = os.path.join(args.path, 'SB{}'.format(args.sbid), 'candidates', 'SB{}.csv'.format(args.sbid))

    # check completion
    peak_path = os.path.join(args.path, 'SB{}'.format(args.sbid), 'candidates')
    peak_num = len([name for name in os.listdir(peak_path) if 'peak_cand.csv' in name])
    print(peak_num)

    if os.path.isfile(fname):
        cand = Table.read(fname)
        try: 
            body = 'Hi, ' + '\n\n' \
                f'SB{args.sbid} finished processing' + '\n' + \
                f'BEAM={peak_num} ' + f'ROW={len(cand)} ' + 'PSR={}'.format(sum(~cand['PSR_name'].mask)) + '\n\n' \
                'The csv should be uploded to google drive!'
        except:
            body = 'Hi, ' + '\n\n' \
                f'SB{args.sbid} finished processing' + '\n' + \
                f'BEAM={peak_num} ' + f'ROW={len(cand)} ' + '\n\n' \
                'The csv should be uploded to google drive!'
        os.system(f'rclone copy {fname} google: -P')
        log.info(body)

    else:
        log.warning("%s doesn't exist", fname)
        sys.exit()


    if args.email is not None:
        # sender_email = input("Enter sender gmail address: ")
        sender_email = os.getenv('SMTP_USER')
        sender_pwd = os.getenv('SMTP_PWD')

        for receiver_email in args.email:
            email_alert(body, sender_email, sender_pwd, receiver_email, args)
            log.info('Yayy!!!')




def email_alert(body, sender_email, sender_pwd, receiver_email, args):

    msg = EmailMessage()
    msg.set_content(body)
    msg['subject'] = f'VASTER SB{args.sbid} finished processing'
    msg['from'] = sender_email
    msg['to'] = receiver_email

    # creates SMTP session
    s = smtplib.SMTP('smtp.gmail.com', 587)
    # start TLS for security
    s.starttls()
    # Authentication
    s.login(sender_email, sender_pwd)
    # sending the mail
    # s.sendmail(sender_email, receiver_email, msgs)
    s.send_message(msg)
    # terminating the session
    s.quit()



if __name__ == '__main__':
    _main()
