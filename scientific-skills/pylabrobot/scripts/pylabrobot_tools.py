#!/usr/bin/env python3
"""
PyLabRobot Tools - Laboratory automation protocol utilities
Liquid handling, plate reading, and equipment control
"""

import argparse
import json
import asyncio
from typing import List, Dict, Optional


def list_supported_devices() -> List[Dict]:
    """List supported laboratory devices."""
    devices = [
        {"type": "Liquid Handler", "models": ["Hamilton STAR", "Hamilton STARlet", "Opentrons OT-2", "Tecan EVO"]},
        {"type": "Plate Reader", "models": ["BMG CLARIOstar"]},
        {"type": "Heater Shaker", "models": ["Hamilton HeaterShaker", "Inheco ThermoShake"]},
        {"type": "Incubator", "models": ["Inheco Incubator", "Thermo Fisher Incubator"]},
        {"type": "Centrifuge", "models": ["Agilent VSpin"]},
        {"type": "Pump", "models": ["Cole Parmer Masterflex"]},
    ]
    return devices


def list_labware_types() -> List[Dict]:
    """List common labware types."""
    labware = [
        {"name": "96-well plate", "wells": 96, "volume": "200-400 uL"},
        {"name": "384-well plate", "wells": 384, "volume": "20-100 uL"},
        {"name": "PCR plate", "wells": 96, "volume": "0.2 mL"},
        {"name": "Deep well plate", "wells": 96, "volume": "1-2 mL"},
        {"name": "Tip rack 1000uL", "tips": 96, "volume": "1000 uL"},
        {"name": "Tip rack 200uL", "tips": 96, "volume": "200 uL"},
        {"name": "Tip rack 10uL", "tips": 96, "volume": "10 uL"},
        {"name": "Reservoir 12-channel", "channels": 12, "volume": "15-20 mL"},
        {"name": "Trough 1-channel", "channels": 1, "volume": "100+ mL"},
    ]
    return labware


def generate_liquid_handling_protocol(output_path: str):
    """Generate a liquid handling protocol template."""
    template = '''#!/usr/bin/env python3
"""
PyLabRobot Liquid Handling Protocol
Example: Sample transfer protocol
"""

import asyncio
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends import ChatterBoxBackend
from pylabrobot.resources import (
    STARLetDeck,
    TIP_CAR_480_A00,
    PLT_CAR_L5AC_A00,
    Cos_96_DW_1mL
)

async def main():
    # Initialize liquid handler with simulation backend
    backend = ChatterBoxBackend()
    lh = LiquidHandler(backend=backend, deck=STARLetDeck())
    await lh.setup()
    
    # Define resources
    tip_rack = TIP_CAR_480_A00(name="tip_rack")
    source_plate = Cos_96_DW_1mL(name="source")
    dest_plate = Cos_96_DW_1mL(name="dest")
    
    # Assign resources to deck
    lh.deck.assign_child_resource(tip_rack, rails=1)
    lh.deck.assign_child_resource(source_plate, rails=10)
    lh.deck.assign_child_resource(dest_plate, rails=15)
    
    # Protocol: Transfer 100 uL from source to destination
    print("Starting protocol...")
    
    # Pick up tips
    await lh.pick_up_tips(tip_rack["A1:H1"])
    
    # Aspirate from source
    await lh.aspirate(source_plate["A1:H1"], vols=100)
    
    # Dispense to destination
    await lh.dispense(dest_plate["A1:H1"], vols=100)
    
    # Drop tips
    await lh.drop_tips()
    
    print("Protocol complete!")

if __name__ == "__main__":
    asyncio.run(main())
'''
    
    with open(output_path, 'w') as f:
        f.write(template)
    
    print(f"Liquid handling protocol saved to {output_path}")


def generate_plate_reading_protocol(output_path: str):
    """Generate a plate reading protocol template."""
    template = '''#!/usr/bin/env python3
"""
PyLabRobot Plate Reading Protocol
Example: Absorbance reading protocol
"""

import asyncio
from pylabrobot.plate_reading import PlateReader
from pylabrobot.plate_reading.clario_star_backend import CLARIOstarBackend

async def main():
    # Initialize plate reader
    reader = PlateReader(name="CLARIOstar", backend=CLARIOstarBackend())
    await reader.setup()
    
    # Set temperature
    print("Setting temperature to 37°C...")
    await reader.set_temperature(37)
    
    # Open plate reader
    await reader.open()
    
    # (Manually or robotically load plate here)
    input("Load plate and press Enter to continue...")
    
    # Close and read
    await reader.close()
    
    # Read absorbance at 450 nm
    print("Reading absorbance...")
    data = await reader.read_absorbance(wavelength=450)
    
    print(f"Results: {data}")
    
    # Save results
    with open("plate_reading_results.json", 'w') as f:
        json.dump(data, f, indent=2)
    
    # Open to remove plate
    await reader.open()
    print("Remove plate")

if __name__ == "__main__":
    asyncio.run(main())
'''
    
    with open(output_path, 'w') as f:
        f.write(template)
    
    print(f"Plate reading protocol saved to {output_path}")


def generate_serial_dilution_protocol(output_path: str):
    """Generate a serial dilution protocol template."""
    template = '''#!/usr/bin/env python3
"""
PyLabRobot Serial Dilution Protocol
Creates 2-fold serial dilutions across a row
"""

import asyncio
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends import ChatterBoxBackend
from pylabrobot.resources import STARLetDeck, TIP_CAR_480_A00, Cos_96_DW_1mL

async def serial_dilution(lh: LiquidHandler, plate, row: str, 
                          dilution_factor: float = 2.0, num_wells: int = 11):
    """
    Perform serial dilution across a row.
    
    Args:
        lh: LiquidHandler instance
        plate: Plate resource
        row: Row letter (A-H)
        dilution_factor: Dilution factor (default 2 for 2-fold)
        num_wells: Number of dilutions to perform
    """
    wells = [f"{row}{i}" for i in range(1, num_wells + 2)]  # A1, A2, ..., A12
    
    # Volume settings
    stock_volume = 100  # uL
    diluent_volume = stock_volume * (dilution_factor - 1)
    transfer_volume = stock_volume
    
    print(f"Performing {dilution_factor}-fold serial dilution on row {row}")
    
    # Add diluent to all wells except first
    for i in range(1, len(wells)):
        await lh.pick_up_tips()
        await lh.aspirate(plate[wells[i]], vols=diluent_volume)
        await lh.dispense(plate[wells[i]], vols=diluent_volume)
        await lh.drop_tips()
    
    # Serial transfer
    for i in range(len(wells) - 1):
        await lh.pick_up_tips()
        await lh.aspirate(plate[wells[i]], vols=transfer_volume)
        await lh.dispense(plate[wells[i + 1]], vols=transfer_volume, 
                         jet=True, blow_out=True)
        await lh.drop_tips()

async def main():
    # Initialize
    backend = ChatterBoxBackend()
    lh = LiquidHandler(backend=backend, deck=STARLetDeck())
    await lh.setup()
    
    # Setup deck
    tip_rack = TIP_CAR_480_A00(name="tips")
    plate = Cos_96_DW_1mL(name="dilution_plate")
    
    lh.deck.assign_child_resource(tip_rack, rails=1)
    lh.deck.assign_child_resource(plate, rails=10)
    
    # Perform serial dilution on row A
    await serial_dilution(lh, plate, "A", dilution_factor=2.0, num_wells=11)
    
    print("Serial dilution complete!")

if __name__ == "__main__":
    asyncio.run(main())
'''
    
    with open(output_path, 'w') as f:
        f.write(template)
    
    print(f"Serial dilution protocol saved to {output_path}")


def calculate_liquid_volumes(protocol_steps: List[Dict]) -> Dict:
    """Calculate total liquid volumes needed for a protocol."""
    volumes = {}
    
    for step in protocol_steps:
        liquid = step.get("liquid", "unknown")
        volume = step.get("volume", 0)
        
        if liquid in volumes:
            volumes[liquid] += volume
        else:
            volumes[liquid] = volume
    
    # Add 20% extra for safety
    volumes = {k: round(v * 1.2, 2) for k, v in volumes.items()}
    
    return volumes


def generate_deck_layout(layout_type: str = "basic") -> Dict:
    """Generate a deck layout configuration."""
    layouts = {
        "basic": {
            "rails_1_8": "Tip carriers",
            "rails_9_14": "Plate carriers",
            "rails_15_22": "Troughs and reagents",
            "rails_23_30": "Waste and additional resources"
        },
        "pcr": {
            "rails_1_4": "Tip carriers (10uL, 200uL)",
            "rails_5_8": "PCR plate carrier",
            "rails_9_12": "Reagent reservoir",
            "rails_13_16": "Master mix plate",
            "rails_17_30": "Sample plates"
        },
        "elisa": {
            "rails_1_4": "Tip carriers",
            "rails_5_8": "Assay plate",
            "rails_9_12": "Wash buffer reservoir",
            "rails_13_16": "Reagent troughs",
            "rails_17_20": "Waste",
            "rails_21_30": "Reader/plate hotel"
        }
    }
    
    return layouts.get(layout_type, layouts["basic"])


def estimate_protocol_time(steps: List[Dict]) -> Dict:
    """Estimate protocol execution time."""
    time_estimates = {
        "pick_up_tips": 5,      # seconds
        "aspirate": 10,          # seconds
        "dispense": 10,          # seconds
        "drop_tips": 5,          # seconds
        "move": 3,               # seconds
        "shake": 60,             # seconds
        "incubate": 300,         # seconds (5 min default)
        "read": 30,              # seconds
    }
    
    total_time = 0
    step_times = []
    
    for step in steps:
        action = step.get("action", "move")
        duration = step.get("duration", time_estimates.get(action, 10))
        step_times.append({"action": action, "estimated_time": duration})
        total_time += duration
    
    return {
        "total_seconds": total_time,
        "total_minutes": round(total_time / 60, 2),
        "steps": step_times
    }


def validate_protocol(protocol_json: str) -> Dict:
    """Validate a protocol JSON file."""
    try:
        with open(protocol_json, 'r') as f:
            protocol = json.load(f)
        
        errors = []
        warnings = []
        
        # Check required fields
        if "name" not in protocol:
            errors.append("Missing protocol name")
        
        if "steps" not in protocol or not protocol["steps"]:
            errors.append("No protocol steps defined")
        
        # Check for common issues
        for i, step in enumerate(protocol.get("steps", [])):
            if "action" not in step:
                errors.append(f"Step {i+1}: Missing action")
            
            if step.get("volume", 0) > 1000:
                warnings.append(f"Step {i+1}: Large volume (>1mL)")
        
        return {
            "valid": len(errors) == 0,
            "errors": errors,
            "warnings": warnings
        }
    except Exception as e:
        return {
            "valid": False,
            "errors": [str(e)],
            "warnings": []
        }


def main():
    parser = argparse.ArgumentParser(description="PyLabRobot Lab Automation Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # List devices command
    subparsers.add_parser("list-devices", help="List supported devices")
    
    # List labware command
    subparsers.add_parser("list-labware", help="List labware types")
    
    # Template commands
    liquid_parser = subparsers.add_parser("template-liquid", help="Generate liquid handling template")
    liquid_parser.add_argument("-o", "--output", default="liquid_handling_protocol.py")
    
    plate_parser = subparsers.add_parser("template-plate", help="Generate plate reading template")
    plate_parser.add_argument("-o", "--output", default="plate_reading_protocol.py")
    
    dilution_parser = subparsers.add_parser("template-dilution", help="Generate serial dilution template")
    dilution_parser.add_argument("-o", "--output", default="serial_dilution_protocol.py")
    
    # Deck layout command
    deck_parser = subparsers.add_parser("deck-layout", help="Show deck layout")
    deck_parser.add_argument("-t", "--type", choices=["basic", "pcr", "elisa"], default="basic")
    
    # Calculate volumes command
    vol_parser = subparsers.add_parser("calculate-volumes", help="Calculate liquid volumes")
    vol_parser.add_argument("protocol", help="Protocol JSON file")
    
    # Estimate time command
    time_parser = subparsers.add_parser("estimate-time", help="Estimate protocol time")
    time_parser.add_argument("protocol", help="Protocol JSON file")
    
    # Validate command
    validate_parser = subparsers.add_parser("validate", help="Validate protocol")
    validate_parser.add_argument("protocol", help="Protocol JSON file")
    
    args = parser.parse_args()
    
    if args.command == "list-devices":
        devices = list_supported_devices()
        print(json.dumps(devices, indent=2))
    
    elif args.command == "list-labware":
        labware = list_labware_types()
        print(json.dumps(labware, indent=2))
    
    elif args.command == "template-liquid":
        generate_liquid_handling_protocol(args.output)
    
    elif args.command == "template-plate":
        generate_plate_reading_protocol(args.output)
    
    elif args.command == "template-dilution":
        generate_serial_dilution_protocol(args.output)
    
    elif args.command == "deck-layout":
        layout = generate_deck_layout(args.type)
        print(json.dumps(layout, indent=2))
    
    elif args.command == "calculate-volumes":
        with open(args.protocol, 'r') as f:
            protocol = json.load(f)
        volumes = calculate_liquid_volumes(protocol.get("steps", []))
        print(json.dumps(volumes, indent=2))
    
    elif args.command == "estimate-time":
        with open(args.protocol, 'r') as f:
            protocol = json.load(f)
        time_estimate = estimate_protocol_time(protocol.get("steps", []))
        print(json.dumps(time_estimate, indent=2))
    
    elif args.command == "validate":
        result = validate_protocol(args.protocol)
        print(json.dumps(result, indent=2))
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
