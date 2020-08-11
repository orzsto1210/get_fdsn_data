from obspy.clients.fdsn import Client
from obspy.clients import iris
from obspy import UTCDateTime
import numpy as np
from geopy import distance

def get_events():
    client = Client("IRIS")
    starttime = UTCDateTime("2020-07-28")
    endtime = UTCDateTime("2020-07-31")
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
    client = Client("IRIS")
    networks = client.get_stations(network="*",station="*",starttime=OT,endtime=OT+120)
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
    client = Client("IRIS")
    bulk = []
    
    for sta in sta_list:
        # P波到時
        client1 = iris.Client()
        result = client1.traveltime(phases=['p'], evloc=(OLat, OLon),
            staloc=[(sta['latitude'], sta['longitude'])], evdepth=ODep)
        # print(result.decode())
        r = result.decode().split()
        p_arrival = OT + float(r[27])
        print(p_arrival)

        b = (sta['network'], sta['station'], "*", "BH?", OT, OT+300)
        bulk.append(b)

    # print(len(bulk))
    
    if bulk:
        data = client.get_waveforms_bulk(bulk, attach_response=True)
        data.remove_response(output="ACC")
        data.resample(100.0)
        print(data)
        # print(len(data))

        # for i in data:
        #     with open(i.get_id()+".txt", "w") as f:
        #         for j in i.data:
        #             f.write(str(j)+"\n")

        # for i in data:
        #     fname=i.get_id()+".png"
        #     i.plot(outfile=fname,format="png")
    

events = get_events()
print(len(events))

for event in events:
    OT, OLon, OLat, ODep, OMag = get_event_info(event)
    sta_list = get_stations(OT, OLon, OLat, ODep, OMag)
    # print(sta_list)
    get_waveforms(sta_list, OT, OLon, OLat, ODep, OMag)
