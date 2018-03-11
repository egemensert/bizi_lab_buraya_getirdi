import numpy as np
import matplotlib.pyplot as plt

to_num = {'i': 1, 'ii': 2, 'iii': 3, 'iv': 4, 'v': 5, 'vi': 6, 'vii': 7}

def plot(l, vil, vil_i, text):
    return None
    plt.plot(vil, l[0][vil_i], 'r.')
    plt.plot([0, vil], [l[0][vil_i], l[0][vil_i]], 'g--', linewidth=0.5)
    plt.plot([vil, vil], [0, l[0][vil_i]], 'g--', linewidth=0.5)
    plt.text(vil -.1, -0.25, text, fontsize=8)
    plt.text(-.3, l[0][vil_i] -.08, "%.2f" % l[0][vil_i], fontsize=8)

def plot_vm(l, vil, text):
    return None
    plt.plot(vil, vil, 'r.')
    plt.plot([0, vil], [vil, vil], 'g--', linewidth=0.5)
    plt.plot([vil, vil], [0, vil], 'g--', linewidth=0.5)
    plt.text(vil -.1, -0.25, text, fontsize=8)
    plt.text(-.3, vil -.08, str(vil), fontsize=8)


class ESpice:
    """ EgemenSPICE for the rescue. """
    T_PAST = 0.
    V_PAST = 0.
    EPS = 1E-6

    def __init__(self, flist):
        self.flist = flist
        self.data = {}
        self._parse()
        self.res = {}

    def calc_data(self, vpath, tpath):
        vp = open(vpath, 'w')
        tp = open(tpath, 'w')
        for part in self.data:
            p = to_num[part]
            voltage_char = []
            for mode in self.data[part]:
                if mode == 'vin':
                    l = self.data[part][mode]
                    vil, vih, vil_i, vih_i = self._get_vil_vih(l)
                    vol, voh, vol_i, voh_i = self._get_vol_voh(l)
                    vm, vm_i = self._get_vm(l)
                    nmh = self._get_nmh(voh, vih)
                    nml = self._get_nml(vil, vol)
                    s = "\t".join(map(str, [p, vil, vih, vol, voh, vm, nmh, nml]))
                    vp.write(s + '\n')
                    plt.title('Q3%s V_{out} vs V_{in} Plot' % part)
                    plt.ylabel('V_{out} (Volts)')
                    plt.xlabel('V_{in} (Volts)')
                    plt.grid()
                    plt.plot(l[1], l[0], 'b-')
                    plot(l, vil, vil_i, 'V_IL')
                    plot(l, vih, vih_i, 'V_IH')
                    plot_vm(l, vm, 'Vm')
                    plt.ylim(0, max(l[1])+0.04)
                    plt.xlim(0, max(l[1]))
                    plt.savefig('plots/q3%s.png' % part)
                    plt.close()
                else:
                    for freq in self.data[part][mode]:
                        f = freq * 1000
                        for cap in self.data[part][mode][freq]:
                            c = cap * 1e-12
                            l = self.data[part][mode][freq][cap]
                            ft, rt, mft, mrt = self._get_times(l, f)
                            mof = self._get_mof(rt, ft)
                            lhp, hlp = self._get_lhp_hlp(mft, mrt, f)
                            spd = self._get_static_power_dissipation()
                            dpd = self._get_dynamic_power_dissipation(f, c, l)
                            mof /= 1e6 # in MHz
                            s = "\t".join(map(str, [p, freq, cap, ft, rt, mof, lhp, hlp, spd, dpd]))
                            tp.write(s + '\n')
    def _parse(self):
        for f in flist:
            fname = f.split('/')[-1]
            if 'txt' not in fname:
                continue
            details = fname[1:-4].split('_')
            part, misc = details[0], details[1]
            if part not in self.data:
                self.data[part] = {
                    'freq':
                        { # freq values
                            1: {1: [], 10: [], 100:[]}, # cap values
                            10: {1: [], 10: [], 100:[]}, # cap values
                            100: {1: [], 10: [], 100:[]}, # cap values
                        },
                    'vin': []}
            if 'khz' in misc:
                mode = 'freq'
                freq = int(misc[:-3])
                logger = self.data[part][mode][freq]
            else:
                mode = 'vin'
                logger = self.data[part][mode]
            with open(f, 'r') as foo:
                address = logger
                cap_index = 0
                for line in foo:
                    code = line.rstrip().split()
                    if code[0] in ['time', 'vin', 'v2']:
                        continue
                    elif code[0] == '0.000000000000000e+000':
                        cap = [1, 10, 100][cap_index]
                        cap_index += 1
                        if mode == 'freq':
                            address = logger[cap]
                    else:
                        address.append(list(map(float, code)))
        for part in self.data:
            for mode in self.data[part]:
                if mode == 'vin':
                    vin, vout = self._wrap(self.data[part][mode])
                    self.data[part][mode] = [vout, vin, self._deriv(vout, vin)]
                    self.data[part][mode] = self._truncate_v(self.data[part][mode])
                    continue
                for freq in self.data[part][mode]:
                    for cap in self.data[part][mode][freq]:
                        t, v = self._wrap(self.data[part][mode][freq][cap])
                        self.data[part][mode][freq][cap] = [v, t, self._deriv(v, t)]
                        self.data[part][mode][freq][cap] = self._truncate_to_period(
                            self.data[part][mode][freq][cap],
                            freq * 1000)

    def _get_static_power_dissipation(self):
        return 0

    def _get_dynamic_power_dissipation(self, freq, cap, l):
        [vo, vi, d] = l
        mv = np.max(vo)
        return freq * cap * mv**2

    def _get_lhp_hlp(self, mft, mrt, freq):
        period = 1. / freq
        hlp = mft
        lhp = mrt - 0.5 * period
        return lhp, hlp

    def _get_mof(self, rise_time, fall_time):
        return 0.5 / max(rise_time, fall_time)

    def _get_times(self, l, freq):
        [v, t, d] = l
        period = 1. / freq
        if max(v) > 6:
            mv = 12.
        else:
            mv = 5.
        fmindex = np.argmin(abs(v[:len(v)//2]- mv * 0.9))
        fmaxdex = np.argmin(abs(v[:len(v)//2]- mv * 0.1))

        rmaxdex  = np.argmin(abs(v[len(v)//2:]- mv * 0.9))
        rmindex  = np.argmin(abs(v[len(v)//2:]- mv * 0.1))
        fall_time = t[fmaxdex] - t[fmindex]
        mid_fall_time = 0.5 * (t[fmaxdex] + t[fmindex])
        rise_time = t[rmaxdex + len(v)//2] - t[rmindex + len(v)//2]
        mid_rise_time = 0.5 * (t[rmaxdex + len(v)//2] + t[rmindex + len(v)//2])
        return fall_time, rise_time, mid_fall_time, mid_rise_time

    def _get_nmh(self, voh, vih):
        return vih - voh

    def _get_nml(self, vol, vil):
        return vol - vil

    def _get_vol_voh(self, l):
        [vo, vi, d] = l
        m = min(vo)
        mi = np.argmin(vo)
        mai = np.argmax(vo)
        if m < self.EPS:
            m = 0
        return m, max(vo), mi, mai

    def _get_vm(self, l):
        [vo, vi, d] = l
        index = np.argmin(abs(vo - vi))
        return vi[index], index

    def _get_vil_vih(self, l):
        [vo, vi, d] = l
        foo = 1 + d
        indexes, = np.where(foo < self.EPS)
        low, high = self._separate_indexes(indexes, d)
        vil = vi[low]
        vih = vi[high]
        return vil, vih, low, high

    def _separate_indexes(self, indexes, d):
        low = []
        high = []
        for i, index in enumerate(indexes):
            if index < len(d) / 2:
                low.append(index)
            else:
                high = indexes[i:]
                break
        low_i = int(np.median(low))
        high_i = int(np.median(high))
        return (low_i, high_i)

    def _wrap(self, x):
        a, b = tuple(zip(*x))
        return np.array(a), np.array(b)

    def _deriv(self, y, x):
        # dy / dx
        y = [y[0]] + y
        x = [x[0]] + x

        dy = self._sliding_difference(y)
        dx = self._sliding_difference(x)

        return dy / (dx + self.EPS)

    def _sliding_difference(self, x):
        foo = []
        for i in range(len(x)-1):
            foo.append(x[i+1] - x[i])
        return np.array(foo)

    def _truncate_to_period(self, l, freq):
        [v, t, d] = l
        period = 1. / freq
        for i in range(len(t)):
            if t[i] > period:
                break
        return [v[:i], t[:i], d[:i]]

    def _truncate_v(self, l):
        [vo, vi, q] = l
        d = max(vi)
        for i, e in enumerate(vi):
            if e == d:
                break
        return [vo[:i+1], vi[:i+1], q[:i+1]]

if __name__ == '__main__':
    import glob
    flist = glob.glob('./*')
    p = ESpice(flist)
    p.calc_data('voltage_char.tsv', 'time_power_char.tsv')
