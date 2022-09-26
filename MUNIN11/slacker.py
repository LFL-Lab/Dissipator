import os
import sys
sys.path.append(r"C:\Users\lfl\measurements_test")

from Slack_wrapper.Slack import Slack_er

slack_channel = "hpc-notifications"
sacha_id = "W018J5SNCCT"

def sendslack(slack_id = f"<@{sacha_id}> ", message = "Measurement finished."):
    Slacker = Slack_er()
    Slacker.send_message(slack_channel, slack_id + message)