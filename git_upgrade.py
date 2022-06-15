#!/usr/bin/env python3
import os
import random
import subprocess
from datetime import datetime

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)


def read_year():
    with open('number.txt', 'r') as f:
        return int(f.read().strip())


def write_year(year):
    with open('number.txt', 'w') as f:
        f.write(str(year))


def git_commit():
    subprocess.run(['git', 'add', 'number.txt'])
    date = datetime.now().strftime('%Y-%m-%d')
    commit_message = f"Year update: {date}"
    subprocess.run(['git', 'commit', '-m', commit_message])


def git_push():
    result = subprocess.run(['git', 'push'], capture_output=True, text=True)
    if result.returncode == 0:
        print("Changes pushed to GitHub successfully.")
    else:
        print("Error pushing to GitHub:")
        print(result.stderr)


def update_cron_with_random_time(year):
    # Commit frequency based on the year
    if year == 2022:
        frequency = 2  # Twice a week
    elif year == 2023:
        frequency = 3  # Three times a week

    # Generate random cron entries for the week
    cron_entries = []
    for _ in range(frequency):
        random_hour = random.randint(0, 23)
        random_minute = random.randint(0, 59)
        random_day_of_week = random.randint(0, 6)  # 0 is Sunday, 6 is Saturday
        cron_command = f"{random_minute} {random_hour} * * {random_day_of_week} cd {script_dir} && python3 {os.path.join(script_dir, 'update_number.py')}\n"
        cron_entries.append(cron_command)

    # Apply to crontab
    cron_file = "/tmp/current_cron"
    os.system(f"crontab -l > {cron_file} 2>/dev/null || true")
    with open(cron_file, "r") as file:
        lines = file.readlines()
    with open(cron_file, "w") as file:
        for line in lines:
            if "update_number.py" not in line:
                file.write(line)
        for entry in cron_entries:
            file.write(entry)
    os.system(f"crontab {cron_file}")
    os.remove(cron_file)
    print(f"Cron job updated for {frequency} times a week in {year}.")


def main():
    try:
        current_year = read_year()
        new_year = 2023 if current_year == 2022 else 2022
        write_year(new_year)
        git_commit()
        git_push()
        update_cron_with_random_time(new_year)
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)


if __name__ == "__main__":
    main()

