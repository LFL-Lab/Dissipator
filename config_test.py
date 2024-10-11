import numpy as np 
from scipy.signal.windows import gaussian
from Utilities import IQ_imbalance
from typing import Dict

class Configuration:
    '''Class to create a configuration file based on the Qsetup parameters'''

    def __init__(self, qsetup):
        self._pars = qsetup.pars
        self._elements = qsetup.pars['elements']
        self._operations = qsetup.pars['operations']
        self._pulses = ['const_pulse', 'X180_pulse', 'Y180_pulse', 'X90_pulse', 'Y90_pulse', 'gaussian_pulse', 'readout_pulse', 'arb_pulse']
        self._controller = {}
        self._waveforms = ['zero_wf', 'const_wf', 'const_wf_rr', 'gaussian_wf', 'readout_wf', 'X180_wf_I', 'X180_wf_Q', 'X90_wf_I', 'X90_wf_Q', 'arb_wfm']
        self._integration_weights = [ 'cos', 'sin', 'minus_sin']
        # 'cos_phi', 'sin_phi', 'minus_sin_phi',
        self._config = {}
        self.init_config()

    def init_config(self):
        self._config = {'version': 1, 'controllers': {self._pars['controller']: {}}}
        self.make_config()

    def add_controller(self):
        controller = dict(type='opx1', analog_outputs=self.make_analog_outputs(),
                          digital_outputs={},
                          analog_inputs=self.make_analog_inputs())
            
        return controller
    
    def make_analog_outputs(self, nChannels:int = 2) -> Dict[int, Dict[str, int]]:
        analog_outputs = {}
        channel_idx = 1
        for element in self._pars['elements']:
            for i in range(nChannels):
                analog_outputs[channel_idx] = {"offset": self._pars[f'{element}_mixer_offsets'][i]}
                channel_idx += 1

        return analog_outputs
    
    def make_analog_inputs(self) -> Dict[int, Dict[str, int]]:
        analog_inputs = {}
        analog_inputs[1] = dict(offset=self._pars['analog_input_offsets'][0],gain_db=self._pars['analog_input_gain'])
        analog_inputs[2] = dict(offset=self._pars['analog_input_offsets'][1],gain_db=self._pars['analog_input_gain'])

        return analog_inputs

    def make_element_dict(self):
        element_dict = {}
        for element in self._pars['elements']:
            if element == 'rr':
                element_dict[element] = self.add_res_element(element)
            else:
                element_dict[element] = self.add_qubit_element(element)
        return element_dict
    
    def add_res_element(self, res_name: str = 'rr', operations: list = []):
        '''
        Generic readout element
        '''
        element =  dict(mixInputs = dict(I=self.make_tuple('Iout',element=res_name),
                         Q=self.make_tuple('Qout',element=res_name),
                         lo_frequency=self._pars[f'{res_name}_LO'],
                         mixer=res_name),                
                intermediate_frequency=self._pars[f'{res_name}_IF'],
                outputs= dict(out1=self.make_tuple('Iin'), out2=self.make_tuple('Qin'),),
                time_of_flight=self._pars['tof'],
                smearing= self._pars['smearing'],
                operations = self.make_operation_dict(res_name),  
                )
        
        return element
    
    def add_qubit_element(self, qb_name: str = 'qubit', operations: list = []):
        '''
        Generic qubit element
        '''

        element = dict(mixInputs = dict(I=self.make_tuple(channel='Iout', element=qb_name),
                                     Q=self.make_tuple(channel='Qout', element=qb_name),
                                     lo_frequency=self._pars[f'{qb_name}_LO'],
                                     mixer=qb_name),                
                    intermediate_frequency=self._pars[f'{qb_name}_IF'],
                    digitalInputs={},
                    operations = self.make_operation_dict(qb_name),)   

        return element 
    
    def make_pulse_dict(self):
        '''
        Create a dictionary of pulses.
        '''
        pulse_dict = {}
        for pulse in self._pulses:
            pulse_dict[pulse] = self.add_pulse(pulse_name=pulse)
        return pulse_dict
    
    def get_pulse_length(self, pulse_name):
        pulse_lengths = {
            'const_pulse': 100,
            'gaussian_pulse': self._pars['gauss_len'],
            'X180_pulse': self._pars['X180_len'],
            'X90_pulse': self._pars['X90_len'],
            'Y180_pulse': self._pars['X180_len'],
            'Y90_pulse': self._pars['X90_len'],
            'arb_pulse': self._pars['arb_op_len'],
            'displace_pulse': self._pars['displace_len'],
        }
        return pulse_lengths.get(pulse_name)

    def get_waveform_names(self, pulse_name):
        waveform_mappings = {
            'const_pulse': dict(I='const_wf', Q='zero_wf'),
            'gaussian_pulse': dict(I='gaussian_wf', Q='zero_wf'),
            'X180_pulse': dict(I='X180_wf_I', Q='X180_wf_Q'),
            'X90_pulse': dict(I='X90_wf_I', Q='X90_wf_Q'),
            'Y180_pulse': dict(I='X180_wf_Q', Q='X180_wf_I'),
            'Y90_pulse': dict(I='X90_wf_Q', Q='X90_wf_I'),
            'arb_pulse': dict(I='arb_wfm', Q='zero_wf'),
            'displace_pulse': dict(I='displace_wf', Q='displace_wf'),
        }
        return waveform_mappings.get(pulse_name)
        
    def add_pulse(self, pulse_name: str = 'const_pulse'):

        if pulse_name not in self._pulses:
            raise ValueError(f'Invalid pulse name! Pulse {pulse_name} not found in list of available pulses')
        else:
            if pulse_name == 'readout_pulse':
                operation_type = 'measurement'
                pulse_length = self._pars['readout_length']
                waveform_names = dict(I = 'readout_wf', Q = 'zero_wf')
            else:
                operation_type = 'control'
                pulse_length = self.get_pulse_length(pulse_name)
                waveform_names = self.get_waveform_names(pulse_name)

        # print(pulse_name,pulse_length)
        assert isinstance(pulse_length, int), 'Pulse length must be an integer'
        assert pulse_length >= 16, 'Pulse length must be greater than 16 ns'

        pulse = dict(operation = operation_type,
                    length = pulse_length,
                    waveforms = waveform_names)
        
        if pulse_name == 'readout_pulse':
            pulse['integration_weights'] = {}
            for weight in self._integration_weights:
                pulse['integration_weights'][weight] = weight
            pulse['digital_marker'] = 'ON'
        else:
            pass
            
        return pulse
                   
    def make_operation_dict(self,element: str):
        operation_dict = {} #initialize
        # print(self._operations[element])
        for operation in self._operations[element]:
            operation_dict[operation] = self.select_pulse_name(element, operation_name=operation)

        return operation_dict

    def select_pulse_name(self, element: str = 'qubit', operation_name:str = 'const'):
    
        if operation_name not in self._pars['operations'][element]:
            raise ValueError(f'Operation {operation_name} not found in operations for {element}')
        else:
            operation_pulse = operation_name + '_pulse'
            # if operation_name == 'const':
            #     operation_pulse = 'const_pulse'
            # elif operation_name == 'gauss':
            #     operation_pulse = 'gaussian_pulse'
            # elif operation_name == 'gauss_4ns':
            #     operation_pulse = 'gaussian_4ns_pulse'
            # elif operation_name == 'X180':
            #     operation_pulse = 'X180_pulse'
            # elif operation_name == 'X90':
            #     operation_pulse = 'X90_pulse'
            # elif operation_name == 'Y180':
            #     operation_pulse = 'Y180_pulse'
            # elif operation_name == 'Y90':
            #     operation_pulse = 'Y90_pulse'
            # elif operation_name == 'readout':
            #     operation_pulse = 'readout_pulse'
            # elif operation_name == 'arb_op':
            #     operation_pulse = 'arb_pulse'    

        return operation_pulse
           
    def make_waveform_dict(self):
        waveform_dict = {}
        for waveform_name in self._waveforms:
            waveform_dict[waveform_name] = self.add_waveform(waveform_name)
        return waveform_dict
    
    def add_waveform(self, wfm_name: str = ''):

        waveform_mapping = {
            "zero_wf": ('constant', 0.0),
            "const_wf": ('constant', self._pars['amp_q']),
            "const_wf_rr": ('constant', self._pars['amp_r']),
            "gaussian_wf": ('arbitrary', [float(arg) for arg in self._pars['gauss_amp'] * gaussian(self._pars['gauss_len'], self._pars['gauss_len']/5)]),
            "readout_wf": ('constant', self._pars['amp_r']),
            "X180_wf_I": ('arbitrary', [float(arg) for arg in self._pars['X180_amp'] * gaussian(self._pars['X180_len'], self._pars['X180_len']/5)]),
            "X180_wf_Q": ('constant', 0.0),
            "X90_wf_I": ('arbitrary', [float(arg) for arg in self._pars['X90_amp'] * gaussian(self._pars['X90_len'], self._pars['X90_len']/5)]),
            "X90_wf_Q": ('constant', 0.0),
            "displace_wf": {"type": "arbitrary", "samples": self._pars['displace_amp'] * gaussian(self._pars['displace_len'], self._pars['displace_sigma'])},
            "arb_wfm": ('arbitrary', [0.2]*self._pars['arb_op_len'])
        }

        if wfm_name not in self._waveforms:
            raise ValueError(f'Invalid waveform name! Waveform {wfm_name} not found in list of available waveforms')
        else:
            waveform_type, waveform_data = waveform_mapping.get(wfm_name)

        if waveform_type == "constant":
            waveform = dict(type=waveform_type, sample=waveform_data)
        elif waveform_type == "arbitrary":
            waveform = dict(type=waveform_type, samples=waveform_data)

        return waveform

        return waveform
    
    def make_integration_weights_dict(self):
            integration_weights_dict = {}
            for weight_name in self._integration_weights:
                integration_weights_dict[weight_name] = self.add_integration_weights(weight_name)
            
            return integration_weights_dict

    def add_integration_weights(self, weight_name):

        if weight_name not in self._integration_weights:
            raise ValueError(f'Invalid integration weight name! Weight {weight_name} not found in list of available weights')
        else:
            weight_mappings = {
                "cos_phi": [(np.cos(self._pars['IQ_rotation']), self._pars['readout_length'])],
                "sin_phi": [(np.sin(self._pars['IQ_rotation']), self._pars['readout_length'])],
                "minus_sin_phi": [(-np.sin(self._pars['IQ_rotation']), self._pars['readout_length'])],
                "cos": [1.0] * self._pars['readout_length'],
                "sin": [0.0] * self._pars['readout_length'],
                "minus_sin": [0.0] * self._pars['readout_length']
            }

            cosine = weight_mappings.get(weight_name, [])
            sine = weight_mappings.get(weight_name, [])

        integration_weights = dict(cosine=cosine, sine=sine)
        
        return integration_weights
    
    def make_mixers_dict(self):
        mixers = {}
        for element in self._pars['elements']:
            mixers[element] = self.add_mixer(element)
        return mixers
    
    def add_mixer(self, element):
        mixer = [dict(intermediate_frequency = self._pars[f'{element}_IF'],
                     lo_frequency = self._pars[f'{element}_LO'],
                     correction = IQ_imbalance(*self._pars[f'{element}_mixer_imbalance']))]
        return mixer
    
    def make_tuple(self, channel:str = 'Iout', element:str = 'res'):
            try:
                chan = (f'{self._pars["controller"]}', self._pars[channel][element])
            except TypeError:
                chan = (f'{self._pars["controller"]}', self._pars[channel])

            return chan
    
    def make_digital_waveform(self):
        dig_wfm = {}
        dig_wfm['ON'] = dict(samples=[tuple((1,0))])
        return dig_wfm
    
    def make_config(self):

        self._config['controllers'][self._pars['controller']] = self.add_controller()
        self._config.update(elements=self.make_element_dict())
        self._config.update(pulses=self.make_pulse_dict())
        self._config.update(waveforms=self.make_waveform_dict())
        self._config.update(digital_waveforms=self.make_digital_waveform())
        self._config.update(integration_weights=self.make_integration_weights_dict())
        self._config.update(mixers=self.make_mixers_dict())

        return self._config
        # print(self._config)
        
    def update_configuration(self, new_pars):
        '''
        Update the configuration class parameters and create a new configuration file with the updated parameters.
        '''
        self._pars = new_pars
        config = self.make_config()

        return config