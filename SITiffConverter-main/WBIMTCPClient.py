import asyncio
from Notebooks.Test_class import Test
from SITiffParser import SITiffParser
from concurrent.futures import ProcessPoolExecutor
import json
import shlex
import subprocess
import sys
from asyncio import Semaphore
# To do: add logger 
class WBIMTCPClient:
    __COMPONENT_NAME = "WBIMTCPClient"
    __DELIMITER = "::"
    __PYTHON_PROCESSING_DELAY_S = 10
    __NUM_SYNC_WORKER = 1

    def __init__(self, host="localhost", port=4000):
        self.host = host
        self.port = port
        self.tasks = []
        self.sync_semaphore = Semaphore(self.__NUM_SYNC_WORKER)

    async def run(self):
        self.reader, self.writer = await asyncio.open_connection(self.host, self.port)
        await self.tcp_client()

    async def write(self, data):
        print(data)
        # out_str = self.__DELIMITER.join((self.__COMPONENT_NAME, f"{data}\n"))
        out_str = f"{self.__COMPONENT_NAME}{self.__DELIMITER}{data}\n"
        self.writer.write(out_str.encode())
        await self.writer.drain()

    async def write_info(self, info, info_data=None):
        if info_data is None: 
            info_str = self.__DELIMITER.join(("INFO", f"{info}")) 
        else:
            info_str = self.__DELIMITER.join(("INFO", f"{info}", f"{info_data}"))
        await self.write(info_str)

    async def write_event(self, event_name, event_data=None):
        if event_data is None:
            env_str = self.__DELIMITER.join(("EVENT", f"{event_name}"))
        else:
            env_str = self.__DELIMITER.join(("EVENT", f"{event_name}", f"{event_data}"))
        # evn_str = f"EVENT{self.__DELIMITER}{event_name}"
        await self.write(env_str)

    async def write_data(self, data_name, data=None):
        if data is None:
            data_str = self.__DELIMITER.join(("DATA", f"{data_name}"))
        else:
            data_str = self.__DELIMITER.join(("DATA", f"{data_name}", f"{data}"))
        await self.write(data_str)

    async def tcp_client(self):
        while True:
            # To be modifivied to handle multi-line command
            data = await self.reader.readline()
            if not data:
                break
            data = data.decode().strip()
            print(data)
            command = data.split(self.__DELIMITER)
            if command[0] == "WBIMTCPServer":
            # assert , "Unexpected command head"
            # Deal with different command here
                if command[1] == "WBIMTCPClient OFF":
                    await self.shut_down_client()
                elif command[1] == "Test":
                    # Create a new process to perform the expensive computation
                    task = asyncio.create_task(
                        self.handle_expensive_computation(Test ,1))
                    task.add_done_callback(self._remove_task)
                    self.tasks.append(task)
                elif command[1] == "ProcessTile":
                    task = asyncio.create_task(self.process_tile(command[2]))
                    task.add_done_callback(self._remove_task)
                    self.tasks.append(task)
                elif command[1] == "SyncTile":
                    msg = json.loads(command[2])
                    print(f"Move folder {msg['source_folder']} to target folder {msg['target_folder']}")
                    await self.sync_folder_to_remote(msg['source_folder'], msg['target_folder'])
                    pass
                elif command[1] == "Setting":
                    await self.handling_settings(command[2], command[3])
                else:
                    # Just send the command back
                    await self.write(f"Unrecognized command:{data}")
            else:
                await self.write(f"Unrecognized command:{data}")

        await self.write("Closing the connection")
        self.writer.close()
        await self.writer.wait_closed()

    def _remove_task(self, task):
        # Check the task state here?
        self.tasks.remove(task)

    async def handle_expensive_computation(self, class_handle, args=(), kwargs=None):
        with ProcessPoolExecutor(max_workers=1) as executor:
            # To be modified here
            # Cannot be simplified to one line. If `process` is called in `__init__` of ExpansiveComputation,
            # the computation will start immediately and blocks the event loop.
            # expensive_computation =
            expensive_computation = ExpensiveComputation(class_handle, args, kwargs)
            result = await asyncio.get_event_loop().run_in_executor(executor, expensive_computation.process)

        # Process the result and send it back to the server
        # To add more handling of the result here
        await self.write(f"Result: {result}")

    async def process_tile(self, tile_info_fp):
        await asyncio.sleep(self.__PYTHON_PROCESSING_DELAY_S)
        print("Finish sleeping")
        with ProcessPoolExecutor(max_workers=1) as executor:
            process_hdl = ExpensiveComputation(SITiffParser, (tile_info_fp, ))
            result = await asyncio.get_event_loop().run_in_executor(executor, process_hdl.process)

        rs_msg = {'doneQ' : result.doneQ}
        rs_msg['tile_info_fp'] = tile_info_fp
        if result.doneQ:
            rs_msg['num_row_shift'] = int(result.num_row_shift) 

        await self.write_event("ProcessTile", json.dumps(rs_msg))

    async def sync_folder_to_remote(self, local_folder, remote_folder):
        # config_file_path = "/"
        async with self.sync_semaphore:
            # Create folder on the remove server:
            # mkdir_command = f"ssh {shlex.quote(remote_folder.split(':')[0])} mkdir -p {shlex.quote(remote_folder.split(':')[1])}"
            # mkdir_args = shlex.split(mkdir_command)

            # mkdir_process = await asyncio.create_subprocess_exec(
            #     *mkdir_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            # )
            # stdout, stderr = await mkdir_process.communicate()
            # if mkdir_process.returncode != 0:
            #         print(f"Failed to create remote directory. Return code: {mkdir_process.returncode}")
            #         print(f"Error: {stderr.decode().strip()}")
            #         return mkdir_process.returncode
            
            # Construct the rsync command
            if (sys.platform == 'win32'):
                wsl_convert_cmd = f"wsl wslpath -a {shlex.quote(local_folder)}"
                local_wls_folder = subprocess.check_output(wsl_convert_cmd, shell=True, text=True).strip()

            command = f"rsync -ra --rsync-path='mkdir -p {shlex.quote(remote_folder.split(':')[1])} && rsync' {shlex.quote(local_wls_folder)} {shlex.quote(remote_folder)}"
            # command = f"rsync -ra {shlex.quote(local_wls_folder)} {shlex.quote(remote_folder)}"
            # command = f"rsync -ra {shlex.quote(local_folder)} {shlex.quote(remote_user)}@{shlex.quote(remote_host)}:{shlex.quote(remote_folder)}"
            # Split the command into arguments for create_subprocess_exec
            args = shlex.split(command)

            if (sys.platform == 'win32'):
                args = ['wsl'] + args
            elif (sys.platform == 'linux'): 
                pass

            # Create and run the subprocess
            process = await asyncio.create_subprocess_exec(
                *args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )

            # Wait for the subprocess to finish and capture the output
            stdout, stderr = await process.communicate()

            # Check if the subprocess completed successfully
            msg = {'return_code' : process.returncode, 'error' : stderr.decode().strip(), \
                   'stdout' : stdout.decode().strip(), 'local_folder' : local_folder, 'command' : command}
            
            await self.write_event("SyncFolder", json.dumps(msg))
            
            return process.returncode      

    async def handling_settings(self, ):
        pass

    async def shut_down_client(self):
        await self.write_info("Waiting for all computational processers to finish...")
        await asyncio.gather(*self.tasks)
        await self.write_info("Closing WBIM TCP client...")
        self.writer.close()
        await self.writer.wait_closed()
        exit()
        
class ExpensiveComputation:
    def __init__(self, class_handle, args=(), kwargs=None):
        self.class_handle = class_handle
        self.args = args
        self.kwargs = kwargs if kwargs is not None else {}
        # self.process()

    def process(self):
        # Perform computationally expensive task
        result = self.class_handle(*self.args, **self.kwargs)
        # result.run()
        return result


if __name__ == "__main__":
    host = "localhost"
    port = 4000
    # This cannot be simplified to a single line. Error: coroutine required in asyncio.run()
    client = WBIMTCPClient(host, port)
    asyncio.run(client.run())
