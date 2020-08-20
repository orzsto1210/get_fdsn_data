from obspy.clients.fdsn import Client
from obspy.clients import iris
from obspy import UTCDateTime, read, Stream, read_events
import numpy as np
from geopy import distance
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

client = Client("IRIS")

def get_events():   
    starttime = UTCDateTime("2020-07-01T00:00:00")
    endtime = UTCDateTime("2020-07-31T23:59:59")
    events = client.get_events(starttime=starttime, endtime=endtime,
                            minmagnitude=6)
    return events

def get_event_info(event):
    OT = event.origins[0].time
    OLon = event.origins[0].longitude
    OLat = event.origins[0].latitude
    ODep = event.origins[0].depth / 1000 # km
    OMag = event.magnitudes[0].mag

    return OT, OLon, OLat, ODep, OMag

def get_networks(OT):
    networks = client.get_stations(network="*",station="*",starttime=OT,endtime=OT+300)
    return networks

def get_station_info(sta):
    return sta.code, sta.latitude, sta.longitude

def get_stations(OT, OLon, OLat, ODep, OMag):
    sta_list = []
    networks = get_networks(OT)
    for net in networks:
        # print("network:", net.code)
        for sta in net:
            sta_name, sta_lat, sta_lon = get_station_info(sta)
            # print("station:", sta_name)
            # print("sta_latitude:", sta_lat)
            # print("sta_longitude:", sta_lon)

            # 計算震央與測站距離有無超過 100 km
            dist = distance.distance((OLat, OLon), (sta_lat, sta_lon)).km
            d = {"network": net.code,
                "station": sta_name,
                "latitude": sta_lat,
                "longitude": sta_lon}
            if dist <= 100:
                sta_list.append(d)
    return sta_list

def get_waveforms(sta_list, OT, OLon, OLat, ODep, OMag):
    bulk = []
    for sta in sta_list:
        b = (sta['network'], sta['station'], "*", "BH?", OT, OT+300)
        bulk.append(b)

    # print(len(bulk))
    try:
        st = client.get_waveforms_bulk(bulk, attach_response=True, minimumlength=290)
        st.remove_response(output="ACC")
        st.resample(100.0)
        # print(st)
        my_st = Stream()
        for i, tr in enumerate(st):
            # print(tr)
            my_st += tr.copy()
            if i % 3 == 2:
                outname = str(tr.stats.starttime)[:19] + "_" + tr.get_id()[:-1] + ".MSEED"
                print(outname)
                my_st.write("mseed/"+outname, format="MSEED")
                my_st = Stream()

        # for tr in st:
        #     with open(tr.get_id()+".txt", "w") as f:
        #         for j in tr.data:
        #             f.write(str(j)+"\n")

        # for tr in st:
        #     fname=tr.get_id()+".png"
        #     tr.plot(outfile=fname,format="png")
    except Exception as e:
        print(e)
    
def get_event_endtime(my_st):
    # merge ENZ
    st = my_st.copy() 
    amp_list = []

    for i in range(len(st[0].data)):
        temp_amp = ((st[0].data[i]**2) + (st[1].data[i]**2) + (st[2].data[i]**2))**(0.5)
        amp_list.append(temp_amp)

    merge_enz = np.asarray(amp_list) 

    # 3軸加速度合成分量連續5秒低於最大PGA的20%
    max_amp = max(merge_enz)
    max_amp_index = np.argmax(merge_enz)
    amp = max_amp * 0.2

    timestamp = 0
    for i in range(max_amp_index, len(merge_enz)-499):
        part = merge_enz[i:i+500]
        count = 0
        for p in part:
            if p < amp:
                count += 1
        if count == 500:
            timestamp = max_amp_index+i
            # print(timestamp)
            break
    # 若找不到結束時間，default=90秒
    if timestamp == 0:
        timestamp = 9000

    return my_st[0].stats.starttime + (timestamp/100)

def get_waveforms_vel(my_st):
    try:
        network = my_st[0].stats.network
        # print(network)
        station = my_st[0].stats.station
        # print(station)
        starttime = my_st[0].stats.starttime
        # print(starttime)
        endtime = get_event_endtime(my_st)
        # print(endtime)
        st = client.get_waveforms(network, station, "*", "BH?", starttime, endtime, attach_response=True)
        st.remove_response(output="VEL")
        st.resample(100.0)
        my_st = Stream()
        for i, tr in enumerate(st):
            # print(tr)
            my_st += tr.copy()
            if i % 3 == 2:
                outname = str(tr.stats.starttime)[:19] + "_" + tr.get_id()[:-1] + ".MSEED"
                print(outname)
                my_st.write("mseed_vel/"+outname, format="MSEED")
                my_st = Stream()
    except Exception as e:
        print(e)

def get_vel_data(my_st, OT, OLon, OLat, ODep, OMag):
    try:
        client1 = iris.Client()

        sta = client.get_stations(network=my_st[0].stats.network, station=my_st[0].stats.station, starttime=OT, endtime=OT+300)
        sta_name, sta_lat, sta_lon = get_station_info(sta[0][0])
        # print(sta_name, sta_lat, sta_lon)
        # P波到時
        result = client1.traveltime(phases=['p', 'P'], evloc=(OLat, OLon),
            staloc=[(sta_lat, sta_lon)], evdepth=ODep)
        # print(result.decode())
        r = result.decode().split()
        sample = 100
        p_arrival = int(float(r[27])*100)

        E = my_st[0].data[p_arrival:p_arrival+sample]
        N = my_st[1].data[p_arrival:p_arrival+sample]
        Z = my_st[2].data[p_arrival:p_arrival+sample]

        pgv, intensity = calc_pgv_intensity(my_st)
        print(pgv, intensity)


    except Exception as e:
        print(e)

def calc_pgv_intensity(my_st):
    st = my_st.copy()

    pgv_list = []
    for i in range(len(st[0].data)):
        tmp_pgv = ((st[0].data[i])**2+(st[1].data[i])**2+(st[2].data[i])**2)**(0.5)
        pgv_list.append(tmp_pgv)
    pgv = max(pgv_list) * 100
    pgv = round(pgv, 2)

    # pgv_z = max(abs(st[0].data))
    # pgv_n = max(abs(st[1].data))
    # pgv_e = max(abs(st[2].data))

    # pgv_z = round(pgv_z, 2)
    # pgv_n = round(pgv_n, 2)
    # pgv_e = round(pgv_e, 2)
    
    if pgv < 0.2:
        intensity = 0
    elif pgv >= 0.2 and pgv < 0.7:
        intensity = 1
    elif pgv >= 0.7 and pgv < 1.9:
        intensity = 2
    elif pgv >= 1.9 and pgv < 5.7:
        intensity = 3
    elif pgv >= 5.7 and pgv < 15:
        intensity = 4
    elif pgv >= 15 and pgv < 30:
        intensity = 5.1
    elif pgv >= 30 and pgv < 50:
        intensity = 5.5
    elif pgv >= 50 and pgv < 80:
        intensity = 6.1
    elif pgv >= 80 and pgv < 140:
        intensity = 6.5
    elif pgv >= 140:
        intensity = 7
    else:
        intensity = -1
        
    return pgv, intensity

def get_waveforms_acc(my_st):
    try:
        network = my_st[0].stats.network
        # print(network)
        station = my_st[0].stats.station
        # print(station)
        starttime = my_st[0].stats.starttime
        # print(starttime)
        endtime = get_event_endtime(my_st)
        # print(endtime)
        st = client.get_waveforms(network, station, "*", "BH?", starttime, endtime, attach_response=True)
        st.remove_response(output="ACC")
        st.resample(100.0)
        my_st = Stream()
        for i, tr in enumerate(st):
            # print(tr)
            my_st += tr.copy()
            if i % 3 == 2:
                outname = str(tr.stats.starttime)[:19] + "_" + tr.get_id()[:-1] + ".MSEED"
                print(outname)
                my_st.write("mseed_acc/"+outname, format="MSEED")
                my_st = Stream()
    except Exception as e:
        print(e)

def get_acc_data(my_st, OT, OLon, OLat, ODep, OMag):
    try:
        client1 = iris.Client()

        sta = client.get_stations(network=my_st[0].stats.network, station=my_st[0].stats.station, starttime=OT, endtime=OT+300)
        sta_name, sta_lat, sta_lon = get_station_info(sta[0][0])
        # print(sta_name, sta_lat, sta_lon)
        # P波到時
        result = client1.traveltime(phases=['p', 'P'], evloc=(OLat, OLon),
            staloc=[(sta_lat, sta_lon)], evdepth=ODep)
        # print(result.decode())
        r = result.decode().split()
        sample = 100
        p_arrival = int(float(r[27])*100)

        E = my_st[0].data[p_arrival:p_arrival+sample]
        N = my_st[1].data[p_arrival:p_arrival+sample]
        Z = my_st[2].data[p_arrival:p_arrival+sample]

        pga, intensity = calc_pga_intensity(my_st)
        print(pga, intensity)

    except Exception as e:
        print(e)

def calc_pga_intensity(my_st):
    st = my_st.copy()
    if max(abs(st[0].data)) == 0:
        pga_z = -1
    else:
        pga_z = max(abs(st[0].data))
    if max(abs(st[1].data)) == 0:
        pga_n = -1
    else:
        pga_n = max(abs(st[1].data))
    if max(abs(st[2].data)) == 0:
        pga_e = -1
    else:
        pga_e = max(abs(st[2].data))
    # print pga_z,pga_n,pga_e
    if pga_z == -1:
        pga = -1
    else:
        pga_list = []
        for i in range(len(st[0].data)):
            tmp_pga = ((st[0].data[i])**2+(st[1].data[i])**2+(st[2].data[i])**2)**(0.5)
            pga_list.append(tmp_pga)
        pga = max(pga_list) *100
        pga = round(pga, 2)

    if pga >= 0.0 and pga < 0.8:
        intensity = 0
    elif pga >= 0.8 and pga < 2.5:
        intensity = 1
    elif pga >= 2.5 and pga < 8.0:
        intensity = 2
    elif pga >= 8 and pga < 25:
        intensity = 3
    elif pga >= 25 and pga < 80:
        intensity = 4
    elif pga >= 80:
        pgv, intensity = calc_pgv_intensity(my_st)
    else:
        intensity = -1
        
    return pga, intensity

# events = get_events()
# events.write("xml/2020-01-01_2020-07-31.xml", format="QUAKEML")
# events = read_events("xml/2020-01-01_2020-07-31.xml")
# print("Total events:", len(events))

#==========Get origin data==========
# for event in tqdm(events):
#     OT, OLon, OLat, ODep, OMag = get_event_info(event)
#     sta_list = get_stations(OT, OLon, OLat, ODep, OMag)
#     # print(sta_list)
#     get_waveforms(sta_list, OT, OLon, OLat, ODep, OMag)

#==========Get Velocity data==========
# files = os.listdir("mseed/")
# files = [x for x in files if x.endswith(".MSEED")]
# print(len(files))

# for i in tqdm(files):
#     my_st = read("mseed/"+i)
#     get_waveforms_vel(my_st)

    # fname=str(my_st[0].stats.starttime)[:19] + "_" + my_st[0].get_id()[:-1] +".png"
    # my_st.plot(outfile="png/"+fname,format="png")

#==========Get Vel data of a few seconds after p arrive==========
# vel_files = os.listdir("mseed_vel/")
# vel_files = [x for x in vel_files if x.endswith(".MSEED")]
# print(len(vel_files))

# for i in tqdm(vel_files):
#     my_st = read("mseed_vel/"+i)
#     # print(my_st[0].stats.starttime, my_st[0].stats.endtime)
#     events = client.get_events(starttime=my_st[0].stats.starttime-1, endtime=my_st[0].stats.endtime)
#     OT, OLon, OLat, ODep, OMag = get_event_info(events[0])
#     get_vel_data(my_st, OT, OLon, OLat, ODep, OMag)

    # fname=str(my_st[0].stats.starttime)[:19] + "_" + my_st[0].get_id()[:-1] +".png"
    # my_st.plot(outfile="png_vel/"+fname,format="png")

#==========Get Acceleration data==========
# files = os.listdir("mseed/")
# files = [x for x in files if x.endswith(".MSEED")]
# print(len(files))

# for i in tqdm(files):
#     my_st = read("mseed/"+i)
#     get_waveforms_acc(my_st)

#==========Get Acc data of a few seconds after p arrive==========
# acc_files = os.listdir("mseed_acc/")
# acc_files = [x for x in acc_files if x.endswith(".MSEED")]
# print(len(acc_files))

# for i in tqdm(acc_files):
#     my_st = read("mseed_acc/"+i)
#     events = client.get_events(starttime=my_st[0].stats.starttime-1, endtime=my_st[0].stats.endtime)
#     OT, OLon, OLat, ODep, OMag = get_event_info(events[0])
#     get_acc_data(my_st, OT, OLon, OLat, ODep, OMag)

    # fname=str(my_st[0].stats.starttime)[:19] + "_" + my_st[0].get_id()[:-1] +".png"
    # my_st.plot(outfile="png_acc/"+fname,format="png")
